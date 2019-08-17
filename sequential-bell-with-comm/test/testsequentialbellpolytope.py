'''
Created on 8 ago. 2019

@author: gsenno
'''
import unittest
from sequentialbellpolytope import SequentialBellPolytope
import qutip as qt
from itertools import product
from sequentialbellscenario import SequentialBellScenario
from behaviour import Behaviour
from qutipauxfunc import createQubitObservable, projectorsForQubitObservable,\
    createMaxEntState


class Test(unittest.TestCase):
    
    def testNumberOfVerticesForTwoAliceAndOneBobDoingChained3Is32768(self):
        poly = SequentialBellPolytope(SequentialBellScenario([[2,2,2],[2,2,2]],[2,2,2]))
        vertices=poly.getListOfVertices()
        self.assertEqual(len(vertices),32768)
    
    def testNumberOfVerticesForTwoAlicesAndOneBobDoingCHSHIs256(self):
        poly = SequentialBellPolytope(SequentialBellScenario([[2,2],[2,2]],[2,2]))
        vertices=poly.getListOfVertices()
        self.assertEqual(len(vertices),256)
    
    def testPRCorrBetweenAlice1andBobAndLocalOtherwiseIsNotInPolytope(self):
        scenario = SequentialBellScenario([[2,2],[2,2]],[2,2])
        expectedCorrelations = {((x1,x2),y,(a1,a2),b):1/2*int((a2==0)&(x1*y==(a1+b)%2))
                                for ((x1,x2),y,(a1,a2),b) in scenario.getTuplesOfEvents()}
        expectedBehaviour = Behaviour(scenario,expectedCorrelations)
        poly = SequentialBellPolytope(scenario)
        self.assertFalse(poly.contains(expectedBehaviour.getProbabilityList()))
        
    def testSingletBetweenAlice1AndBobAndLocalOtherwiseIsNotInPolytope(self):
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
        poly = SequentialBellPolytope(scenario)
        self.assertFalse(poly.contains(expectedBehaviour.getProbabilityList()))
         
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testStrategiesForFirstAliceDonnotDependOnSecondAliceInput']
    unittest.main()