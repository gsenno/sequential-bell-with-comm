
'''
Created on Feb 12, 2019

@author: rravell
'''
import numpy as np
from qutipauxfunc import computeDistributionFromStateAndEffects
import qutip as qt
from itertools import product, islice
import picos as pic
from bellpolytopewithonewaycomm import BellPolytopeWithOneWayCommunication
from bellpolytope import BellPolytope
from bellscenario import BellScenario

def createEffects(eigenvectorsPerInput):
    return list(map(lambda eigenvectors :  
                   list(map(lambda eigvector : eigvector*eigvector.dag(),eigenvectors)),eigenvectorsPerInput))

if __name__ == '__main__':
    
    K=2
    ratios=[0,0.464,0.569,0.636,0.683,0.718]
    outputsAlice = [2 for _ in range(K+1)]
    outputsBob = [2 for _ in range(K+1)]
    
    ratio = ratios[K]
    beta=1/(1+ratio)
    alpha=beta*ratio
    c = [(-1)**k * beta**(k+1/2)/np.sqrt(beta**(2*k+1)+alpha**(2*k+1)) for k in range(K+1)]
    
    states=[c[k]*qt.basis(2, 0)+np.sqrt(1-c[k]**2)*qt.basis(2, 1)
                  for k in range(K+1)]
    
    psi=(1/np.sqrt(alpha**2+beta**2))*(alpha*qt.tensor(qt.basis(2,0),qt.basis(2,0))-beta*qt.tensor(qt.basis(2,1),qt.basis(2,1)))
    aliceEffects=[[state*state.dag(),qt.qeye(2)-state*state.dag()] for state in states]
    bobEffects=aliceEffects
    
    dist = computeDistributionFromStateAndEffects(psi,aliceEffects,bobEffects)
    print(dist)
    print(sum(dist))
    polytope = BellPolytopeWithOneWayCommunication(BellPolytope(BellScenario(outputsAlice,outputsBob)))                
    
    if not polytope.contains(dist):
        print('We found a qdist not reproducible with one bit of comm!')