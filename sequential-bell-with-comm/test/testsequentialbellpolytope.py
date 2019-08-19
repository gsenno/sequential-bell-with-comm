'''
Created on 8 ago. 2019

@author: gsenno
'''
import unittest
from sequentialbellpolytope import SequentialBellPolytope
import qutip as qt
import numpy as np
from itertools import product
from sequentialbellscenario import SequentialBellScenario
from behaviour import Behaviour
from qutipauxfunc import createQubitObservable, projectorsForQubitObservable,\
    createMaxEntState, effectForQubitPovm
from bellpolytopewithonewaycomm import BellPolytopeWithOneWayCommunication


class Test(unittest.TestCase):
    
    def testNumberOfVerticesForTwoAliceAndOneBobDoingChained3Is32768(self):
        poly = SequentialBellPolytope(SequentialBellScenario([[2,2,2],[2,2,2]],[2,2,2]))
        vertices=poly.getListOfVertices()
        self.assertEqual(len(vertices),32768)
        
    def testNumberOfVerticesForTwoAliceAndOneBobDoingChained3Is2304(self):
        poly = SequentialBellPolytope(SequentialBellScenario([[2,2],[2,2,3]],[2,2]))
        vertices=poly.getListOfVertices()
        self.assertEqual(len(vertices),2304)
        
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
        alice1Krauss = list(map(lambda qubitObservable : projectorsForQubitObservable(qubitObservable),alice1Observables))
    
        bobUnBlochVectors = [[-1,-1,0],[-1,1,0]]
        bobObservables=list(map(lambda bloch : createQubitObservable(bloch),bobUnBlochVectors))
        bobEffects = list(map(lambda qubitObservable : projectorsForQubitObservable(qubitObservable),bobObservables))
    
        psi=createMaxEntState(2)
        
        expectedCorrelations = {((x1,x2),y,(a1,a2),b):int((a2==0))*((qt.tensor(alice1Krauss[x1][a1],bobEffects[y][b])*psi*psi.dag()).tr())
                                for ((x1,x2),y,(a1,a2),b) in scenario.getTuplesOfEvents()}
        expectedBehaviour = Behaviour(scenario,expectedCorrelations)
        poly = SequentialBellPolytope(scenario)
        self.assertFalse(poly.contains(expectedBehaviour.getProbabilityList()))
    
    def testMeasurementsOverSeparableStateAreLocal(self):
        alice1outputs=[2,2]
        alice2outputs=[2,3]
        bobOutputs=[2,2]
        scenario = SequentialBellScenario([alice1outputs,alice2outputs],bobOutputs)
       
        epsilon=0
        plus=1/np.sqrt(2)*(qt.basis(2, 0)+qt.basis(2, 1))
        minus=1/np.sqrt(2)*(qt.basis(2, 0)-qt.basis(2, 1))
        Kplus=np.cos(epsilon)*plus*plus.dag()+np.sin(epsilon)*minus*minus.dag()
        Kminus=-np.cos(epsilon)*minus*minus.dag()+np.sin(epsilon)*plus*plus.dag()
        
        alice1Krauss=[projectorsForQubitObservable(createQubitObservable([0,0,1])),[Kplus,Kminus]]
        
        alice2blochVectors=[[0,0,1]]
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
        
        psi=qt.tensor(qt.basis(2,0),qt.basis(2,0))
        rho=psi*psi.dag()
        expectedCorrelations={}
        for x1,x2,y in product(range(len(alice1outputs)),range(len(alice2outputs)),range(len(bobOutputs))):
            for a1,a2,b in product(range(alice1outputs[x1]),range(alice2outputs[x2]),range(bobOutputs[y])):
                postMeasrmntState = qt.tensor(alice1Krauss[x1][a1],qt.qeye(2))*rho*(qt.tensor(alice1Krauss[x1][a1],qt.qeye(2))).dag()
                expectedCorrelations[(x1,x2),y,(a1,a2),b]=(qt.tensor(alice2Effects[x2][a2],bobEffects[y][b])*
                                                               postMeasrmntState).tr().real
        
        expectedBehaviour = Behaviour(scenario,expectedCorrelations)
    
        polytope = BellPolytopeWithOneWayCommunication(SequentialBellPolytope(scenario))                
        
        self.assertTrue(polytope.contains(expectedBehaviour.getProbabilityList()))
        
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testStrategiesForFirstAliceDonnotDependOnSecondAliceInput']
    unittest.main()