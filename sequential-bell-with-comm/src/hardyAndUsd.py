
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
    
    outputsAlice = [2,2]
    outputsBob = [2,2]
    
    qzero=qt.basis(2,0)
    qone=qt.basis(2,1)
    plus=1/np.sqrt(2)*(qzero+qone)
    minus=1/np.sqrt(2)*(-qzero+qone)
    
    psi=2/np.sqrt(3)*(qt.tensor(qzero,qzero)-1/2*qt.tensor(plus,plus))
    
    aliceEffects=[[qzero*qzero.dag(),qone*qone.dag()],[plus*plus.dag(),minus*minus.dag()]]
    bobEffects=[[qzero*qzero.dag(),qone*qone.dag()],[plus*plus.dag(),minus*minus.dag()]]
    
    bobPostMeasurmentStates=[[1/((qt.tensor(effect,qt.qeye(2))*psi*psi.dag()).tr())*(qt.tensor(effect,qt.qeye(2))*psi).ptrace(1)
                               for effect in measurement] 
                               for measurement in aliceEffects]
    usdState1=bobPostMeasurmentStates[0][0]
    usdState2=bobPostMeasurmentStates[1][1]
    overlap=np.sqrt((usdState1*usdState2).tr())
    usd=[(qt.qeye(2)-usdState1)/(1+overlap),
         (qt.qeye(2)-usdState2)/(1+overlap)]
    
    #bobEffects.append([usd[0],usd[1],qt.qeye(2)-usd[0]-usd[1]])
     
    dist = computeDistributionFromStateAndEffects(psi,aliceEffects,bobEffects)
    dist = [p.real for p in dist]
    
    scenario=BellScenario(outputsAlice,outputsBob)
    poly=BellPolytopeWithOneWayCommunication(BellPolytope(scenario))
    
    print(poly.contains(dist))
