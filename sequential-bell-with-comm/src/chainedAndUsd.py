
'''
Created on Feb 12, 2019

@author: rravell
'''
import numpy as np
import qutip as qt
from itertools import product, islice
from bellpolytopewithonewaycomm import BellPolytopeWithOneWayCommunication
from qutipauxfunc import createMaxEntState, projectorsForQubitObservable,\
    createQubitObservable, computeDistributionFromStateAndEffects
from bellscenario import BellScenario
from bellpolytope import BellPolytope
from behaviour import Behaviour

def createEffects(eigenvectorsPerInput):
    return list(map(lambda eigenvectors :  
                   list(map(lambda eigvector : eigvector*eigvector.dag(),eigenvectors)),eigenvectorsPerInput))

if __name__ == '__main__':
    n=4
    
    outputsAliceSequence = [2 for _ in range(n)]
    outputsBob = [2 for _ in range(n)]+[3]
    
    phi = lambda j : (j-1)*np.pi/n
    phiprime= lambda j : (2*j-1)*np.pi/(2*n) 
    
    psi=createMaxEntState(2)
    
    aliceEffects=[projectorsForQubitObservable
                  (createQubitObservable([np.sin(phi(j)),0,np.cos(phi(j))]))
                     for j in range(1,n+1)]
        
    bobEffects=[projectorsForQubitObservable
                    (createQubitObservable([np.sin(phiprime(j)),0,np.cos(phiprime(j))]))
                     for j in range(1,n+1)]
    
    bobPostMeasurmentStates=[[1/((qt.tensor(effect,qt.qeye(2))*psi*psi.dag()).tr())*(qt.tensor(effect,qt.qeye(2))*psi).ptrace(1)
                               for effect in measurement] 
                               for measurement in aliceEffects]
    usdState1=bobPostMeasurmentStates[0][0]
    usdState2=bobPostMeasurmentStates[3][1]
    overlap=np.sqrt((usdState1*usdState2).tr())
    usd=[(qt.qeye(2)-usdState1)/(1+overlap),
         (qt.qeye(2)-usdState2)/(1+overlap)]
    
    bobEffects.append([usd[0],usd[1],qt.qeye(2)-usd[0]-usd[1]])
     
    dist = computeDistributionFromStateAndEffects(psi,aliceEffects,bobEffects)
    dist = [p.real for p in dist]
    
    scenario=BellScenario(outputsAliceSequence,outputsBob)
    poly=BellPolytopeWithOneWayCommunication(BellPolytope(scenario))
    
    print(poly.contains(dist))
