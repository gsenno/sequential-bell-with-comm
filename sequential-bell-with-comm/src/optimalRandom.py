
'''
Created on Feb 12, 2019

@author: gsenno
'''

import numpy as np
import qutip as qt
from itertools import product
from mosek.fusion import *

import picos as pic

from bellpolytopewithonewaycomm import BellPolytopeWithOneWayCommunication
from sequentialbellpolytope import SequentialBellPolytope
from sequentialbellscenario import SequentialBellScenario
from qutipauxfunc import createQubitObservable, projectorsForQubitObservable,\
    createMaxEntState, effectForQubitPovm
from behaviour import Behaviour


if __name__ == '__main__':
    
    alice1outputs=[2,2]
    alice2outputs=[2,2,3]
    bobOutputs=[2,2]
    scenario = SequentialBellScenario([alice1outputs,alice2outputs],bobOutputs)
   
    epsilon=7*np.pi/32
    plus=1/np.sqrt(2)*(qt.basis(2, 0)+qt.basis(2, 1))
    minus=1/np.sqrt(2)*(qt.basis(2, 0)-qt.basis(2, 1))
    Kplus=np.cos(epsilon)*plus*plus.dag()+np.sin(epsilon)*minus*minus.dag()
    Kminus=-np.cos(epsilon)*minus*minus.dag()+np.sin(epsilon)*plus*plus.dag()
    
    alice1Krauss=[projectorsForQubitObservable(createQubitObservable([0,0,1])),[Kplus,Kminus]]
    
    alice2blochVectors=[[0,0,1],[1,0,0]]
    alice2Observables = list(map(lambda bloch : createQubitObservable(bloch),alice2blochVectors))
    alice2Effects = list(map(lambda qubitObservable : projectorsForQubitObservable(qubitObservable),alice2Observables))
    trineAlice2 = [[0,0,1],[np.sin(2*np.pi/3),0,np.cos(2*np.pi/3)],[np.sin(4*np.pi/3),0,np.cos(4*np.pi/3)]]
    paulies=[qt.sigmax(),qt.sigmay(),qt.sigmaz()]
    alice2Effects.append(list(map(lambda bloch : 
                           effectForQubitPovm(1/3, sum([paulies[i]*bloch[i] for i in range(3)])),trineAlice2)))
    
    mu=np.arctan(np.sin(2*epsilon))
    bobUnBlochVectors = [[np.sin(mu),0,np.cos(mu)],[-np.sin(mu),0,np.cos(mu)]]
    bobObservables=list(map(lambda bloch : createQubitObservable(bloch),bobUnBlochVectors))
    bobEffects = list(map(lambda qubitObservable : projectorsForQubitObservable(qubitObservable),bobObservables))
    
    psi=createMaxEntState(2)
    rho=psi*psi.dag()
    
    expectedCorrelations={}
    for x1,x2,y in product(range(len(alice1outputs)),range(len(alice2outputs)),range(len(bobOutputs))):
        for a1,a2,b in product(range(alice1outputs[x1]),range(alice2outputs[x2]),range(bobOutputs[y])):
            postMeasrmntState = qt.tensor(alice1Krauss[x1][a1],qt.qeye(2))*rho*(qt.tensor(alice1Krauss[x1][a1],qt.qeye(2))).dag()
            expectedCorrelations[(x1,x2),y,(a1,a2),b]=(qt.tensor(alice2Effects[x2][a2],bobEffects[y][b])*
                                                           postMeasrmntState).tr().real
    
    expectedBehaviour = Behaviour(scenario,expectedCorrelations)
    dist=expectedBehaviour.getProbabilityList()
    
    polytope = BellPolytopeWithOneWayCommunication(SequentialBellPolytope(scenario))                
    vertices = polytope.getListOfVertices()
    with Model("lo1") as M:
    
        # Create variable 'x' of length 4
        x = M.variable("x", len(vertices), Domain.greaterThan(0.0))
         
        
        # Create constraints
        for prob in range(len(vertices[0])):
            M.constraint('p('+str(prob)+')',Expr.dot(x,list(map(lambda ver : ver[prob],vertices))),Domain.equalsTo(dist[prob]))
             
        M.constraint('norm',Expr.dot(x,np.ones((len(vertices), 1))),Domain.equalsTo(1))
         
         
        # Set the objective function to (c^t * x)
        M.objective("obj", ObjectiveSense.Maximize, 1)
        
        # Solve the problem
        M.solve()
        
        # Get the solution values
        print(M.getProblemStatus(SolutionType.Basic))