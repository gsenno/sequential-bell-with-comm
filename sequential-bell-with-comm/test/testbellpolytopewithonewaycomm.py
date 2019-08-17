
'''
Created on 5 ago. 2019

@author: gsenno
'''
import unittest
from bellpolytopewithonewaycomm import BellPolytopeWithOneWayCommunication
import numpy as np
from bellpolytope import BellPolytope
from bellscenario import BellScenario


class Test(unittest.TestCase):

    def testNumberOfVerticesForCHSHIs64(self):
        vertices = BellPolytopeWithOneWayCommunication(BellPolytope(BellScenario([2,2],[2,2]))).getListOfVertices()
        self.assertEqual(len(vertices), 64)
    
    def testExistsStrategyAchievingAlgMaxOfCHSH(self):
        vertices = BellPolytopeWithOneWayCommunication(BellPolytope(BellScenario([2,2],[2,2]))).getListOfVertices()
        chshFunctional = [1,-1,-1,1,1,-1,-1,1,1,-1,-1,1,-1,1,1,-1]
        self.assertEqual(max([np.dot(vertice.getProbabilityList(),chshFunctional) for vertice in vertices]), 4)
        
    def testNumberOfVerticesInScenarioWith4BinaryInputsForBobAnd2ForAliceIs1408(self):
        vertices = BellPolytopeWithOneWayCommunication(BellPolytope(BellScenario([2,2],[2,2,2,2]))).getListOfVertices()
        self.assertEqual(len(vertices), 1408)
    
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testNumberOfVerticesForCHSHIs']
    unittest.main()