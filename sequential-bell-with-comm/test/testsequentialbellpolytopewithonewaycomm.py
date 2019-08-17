'''
Created on 5 ago. 2019

@author: gsenno
'''
import unittest
from bellpolytopewithonewaycomm import BellPolytopeWithOneWayCommunication
import numpy as np
from sequentialbellscenario import SequentialBellScenario
from behaviour import Behaviour
from sequentialbellpolytope import SequentialBellPolytope
from linopttools import createQubitObservable, projectorsForQubitObservable,\
    createMaxEntState
import qutip as qt
from itertools import product

class Test(unittest.TestCase):
    
    def testPRCorrBetweenAlice1andBobAndLocalOtherwiseIsInPolytope(self):
        scenario = SequentialBellScenario([[2,2],[2,2]],[2,2])
        expectedCorrelations = {((x1,x2),y,(a1,a2),b):1/2*int((a2==0)&(x1*y==(a1+b)%2))
                                for ((x1,x2),y,(a1,a2),b) in scenario.getTuplesOfEvents()}
        expectedBehaviour = Behaviour(scenario,expectedCorrelations)
        poly = BellPolytopeWithOneWayCommunication(SequentialBellPolytope(scenario))
        self.assertTrue(poly.contains(expectedBehaviour.getProbabilityList()))
        
    def testSingletBetweenAlice1AndBobAndLocalOtherwiseIsInPolytope(self):
        scenario = SequentialBellScenario([[2,2],[2,2]],[2,2])
        alice1blochVectors = [[1,0,0],[0,1,0]]
        alice1Observables = list(map(lambda bloch : createQubitObservable(bloch),alice1blochVectors))
        alice1Effects = list(map(lambda qubitObservable : projectorsForQubitObservable(qubitObservable),alice1Observables))
    
        bobUnBlochVectors = [[-1,-1,0],[-1,1,0]]
        bobObservables=list(map(lambda bloch : createQubitObservable(bloch),bobUnBlochVectors))
        bobEffects = list(map(lambda qubitObservable : projectorsForQubitObservable(qubitObservable),bobObservables))
    
        psi=createMaxEntState(2)
        
        expectedCorrelations = {((x1,x2),y,(a1,a2),b):int((a2==0))*((qt.tensor(alice1Effects[x1][a1],bobEffects[y][b])*psi*psi.dag()).tr())
                                for ((x1,x2),y,(a1,a2),b) in scenario.getTuplesOfEvents()}
        expectedBehaviour = Behaviour(scenario,expectedCorrelations)
        poly = BellPolytopeWithOneWayCommunication(SequentialBellPolytope(scenario))
        self.assertTrue(poly.contains(expectedBehaviour.getProbabilityList()))
    
    def testProjectiveMeasurementsAreSimulable(self):
        scenario = SequentialBellScenario([[2,2],[2,2]],[2,2])
        alice1blochVectors = [[1,0,0],[0,1,0]]
        alice1Observables = list(map(lambda bloch : createQubitObservable(bloch),alice1blochVectors))
        alice1Effects = list(map(lambda qubitObservable : projectorsForQubitObservable(qubitObservable),alice1Observables))
        
        phases=[np.random.uniform(-np.pi/2,np.pi/2),np.random.uniform(-np.pi/2,np.pi/2)]
        alice2blochVectors=[[np.sin(theta),0,np.cos(theta)] for theta in phases]
        alice2Observables = list(map(lambda bloch : createQubitObservable(bloch),alice2blochVectors))
        alice2Effects = list(map(lambda qubitObservable : projectorsForQubitObservable(qubitObservable),alice2Observables))
        
        bobUnBlochVectors = [[-1,-1,0],[-1,1,0]]
        bobObservables=list(map(lambda bloch : createQubitObservable(bloch),bobUnBlochVectors))
        bobEffects = list(map(lambda qubitObservable : projectorsForQubitObservable(qubitObservable),bobObservables))
    
        psi=createMaxEntState(2)
        rho=psi*psi.dag()
        expectedCorrelations={}
        for x1 in range(2):
            for a1 in range(2):
                postMeasrmntState = qt.tensor(alice1Effects[x1][a1],qt.qeye(2))*rho*(qt.tensor(alice1Effects[x1][a1],qt.qeye(2))).dag()
                for x2,y,a2,b in product(range(2),repeat=4):
                    expectedCorrelations[(x1,x2),y,(a1,a2),b]=(qt.tensor(alice2Effects[x2][a2],bobEffects[y][b])*
                                                               postMeasrmntState).tr().real
        expectedBehaviour = Behaviour(scenario,expectedCorrelations)
                                                    
        poly = BellPolytopeWithOneWayCommunication(SequentialBellPolytope(scenario))
        self.assertTrue(poly.contains(expectedBehaviour.getProbabilityList()),'phases:'+str(phases[0])+', '+str(phases[1]))
        
          
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testNumberOfVerticesForCHSHIs']
    unittest.main()