#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 12:42:14 2020

@author: danielstilckfranca
"""
import qutip as qt
import numpy as np 
import random as random



def random_bipartite(n,k):
    U=qt.rand_unitary_haar(n*k,dims=[[n,k], [n,k]])
    return U

def random_state(n):
    #generates random product state of n qubits
    rho1=qt.rand_dm(2)
    for k in range(0,n-1):
        rho1=qt.tensor(rho1,qt.rand_dm(2))
    return rho1





def random_state_basis(n):
    #generates random tensor product of states consisting of +/- and 0,1
    
    x=random.randint(0,1)
    y=random.randint(0,1)
    initial=qt.basis(2, 0)
    #Hadamard gate
    H=qt.qip.operations.snot()
    #sigma_x gate
    X=qt.sigmax()
    
    initial=(X**(x))*(H**(y))*initial
    initial=initial*initial.dag()
    state=initial
    for k in range(0,n-1):
        x=random.randint(0,1)
        y=random.randint(0,1)
        initial=qt.basis(2, 0)
        #Hadamard gate
        H=qt.qip.operations.snot()
        #sigma_x gate
        X=qt.sigmax()
        
        initial=(X**(x))*(H**(y))*initial
        initial=initial*initial.dag()
        state=qt.tensor(initial,state)
        
    return state
    
    

def max_dist(elements):
    #computes the maximum absolute value of a list
    maximum=0
    for x in elements:
        for y in elements:
            maximum=max(maximum,np.abs(x-y))
    return maximum

def shallow_circuit(n,D):
    #generates a random circuit of depth D on n qubits arragend on a line
    layers=[]
    unitaries=[]
    unitaries.append(qt.rand_unitary_haar(2))
    
    for j in range(0,int((n-1)/2)):
        unitaries.append(qt.rand_unitary_haar(4,dims=[[2,2], [2,2]]))
    unitaries.append(qt.rand_unitary_haar(2))
    
    layers.append(unitaries)
    
    
    
    for k in range(0,D):
        unitaries=[]
        if k%2 is 1:
            unitaries.append(qt.rand_unitary_haar(2))
            
            
            for j in range(0,int((n-1)/2)):
                unitaries.append(qt.rand_unitary_haar(4,dims=[[2,2], [2,2]]))
            unitaries.append(qt.rand_unitary_haar(2))
            layers.append(unitaries)
            
        else:
            unitaries.append(qt.rand_unitary_haar(4,dims=[[2,2], [2,2]]))
            for j in range(0,int((n-1)/2)):
                unitaries.append(qt.rand_unitary_haar(4,dims=[[2,2], [2,2]]))

        layers.append(unitaries)
    return layers
    
    
    
def apply_channel(rho,U,sigma):
    
    #applies a quantum channel with Stinespring unitary U and initial state sigma to a state rho
    n_qubits=len(rho.dims[0])
    tau=qt.tensor(rho,sigma)
    
    tau=U*tau*U.dag()
    return tau.ptrace(range(0,n_qubits))

def apply_channel_slow(rho,U,sigma,q):
    #applies a quantum channel with Stinespring unitary U and initial state sigma to a state rho
    #with probability q, otherwise just applies the identity
    n_qubits=len(rho.dims[0])
    if (random.random()<q):
        tau=qt.tensor(rho,sigma)

        tau=U*tau*U.dag()
        
        return tau.ptrace(range(0,n_qubits))
    else:
        return rho
    


def random_pauli(k):
    #returns a random pauli of 1 or 2 qubits
    if k is 1:
        j=random.randint(0,2)
        if j is 0:
            return qt.sigmax()
        if j is 1:
            return qt.sigmaz()
        if j is 2:
            return qt.sigmay()
    else:
        return qt.tensor(random_pauli(1),random_pauli(1))




def noisy_circuit(layers,p):
    #generates a noisy circuit where each gate is followed with a random pauli gate with probability 
    #p. Note that this generates depolarizing noise with probability p in expectation.
    V=layers[0]
    if (random.random()<p):
    
        U=V[0]
    else:
        U=V[0]*random_pauli(len(V[0].dims[0]))

    for k in range(1,len(V)):
        if (random.random()<p):
            
            U=qt.tensor(U,V[k])
        else:
            Pauli=random_pauli(len(V[k].dims[0]))
            U=qt.tensor(U,V[k]*random_pauli(len(V[k].dims[0])))
    Uf=U
    for k in range(1,len(layers)):
        V=layers[k]
        U=V[0]
        for k in range(1,len(V)):
            U=qt.tensor(U,V[k])
            
        Uf=Uf*U
    return Uf



#example script to generate an MPS of 50 qubits and check convergence




#n here is the number of quibts in the bath.  We start by specifying the observable, called obs, we want to 
    #measure in the bath. In this case, we will measure sigma_x on the first two qubits of the bath.
#note that the code only works for n odd due to the assignment of the gates.
n=5

obs=qt.sigmax()
obs=qt.tensor(qt.sigmax(),obs)


for k in range(0,n-2):

    obs=qt.tensor(qt.qeye(2),obs)


#we will now set the inital state of the bath to be the 0 state.
initial=qt.basis(2, 0)
initial=initial*initial.dag()
right=initial

for k in range(0,n-1):
    right=qt.tensor(initial,right)
    

#we will now generate the rnadom circuit on the bath.
#depth is the depth of the circuit we are going to generate.
depth=3
U=shallow_circuit(n,depth)

#sigma will be the initial state of the qubits coming in to interact with the bath, which we pick at random.
#note that sigma is fixed throughout the run, i.e. all incoming states are the same.
sigma=qt.rand_dm(2)    

#this will apply the quantum channel defined via U and sigma 50 times to determine the expectation vallue
#of running the whole circuit. 
#we choose p=1, so that there
#is 0 depolarizing probability
p=1
right_noise=right

#convert layers to unitary
V=noisy_circuit(U,p)
for k in range(0,50):

    right_noise=apply_channel(right_noise,V,sigma)


#this gives the expectation value if we ran the whole circuit
final=qt.expect(obs,right_noise)


#now we set the number of samples m we want to generate to estimate convergence
m=10
samples=m



#max_length sets how far back we want to go. The results_length list stores the 
#results for a given length. Results stores the samples from that length
results_length=[]
max_length=15

#we can also reset p to compare the value of the circuit to noise.
p=0.99

for k in range(0,max_length):
    results=[]
    for j in range(0,samples):
            
        rho1=random_state_basis(n)
        for length in range(0,k):
                
            rho1=apply_channel(rho1,noisy_circuit(U,p),sigma)
        results.append(qt.expect(rho1,obs))
    results_length.append(results)


#the results are then stored in results length. we can then evaluate the maximal distance of two
    #points on the list by using the max_dist function.
    


