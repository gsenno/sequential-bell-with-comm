'''
Created on 5 ago. 2019

@author: gsenno
'''
from bellpolytope import BellPolytope
from itertools import product
from behaviour import Behaviour

'''
Polytope for bipartite strategies using shared randomness and 1 bit of communication from Bob to Alice.
'''
class BellPolytopeWithOneWayCommunication(BellPolytope):
    
    def getBellScenario(self):
        return self.underlyingPolytope.getBellScenario()
    
    #see bellpolytope.py for description of params
    def __init__(self, underlyingPolytope):
        self.underlyingPolytope=underlyingPolytope
        
    def getGeneratorForVertices(self):
        #local vertices
        yield from self.underlyingPolytope.getGeneratorForVertices()
        
        aliceStrategiesWithComm = [(a1,a2) 
                                   for a1,a2 in product(self.underlyingPolytope.getAliceStrategies(),repeat=2)
                                   if not a1==a2]
         
        communicationStrgs=[format(i,'0'+str(self.underlyingPolytope.numberOfInputsBob())+'b') 
                            for i in range(1,2**(self.underlyingPolytope.numberOfInputsBob()-1))]
                                                    
        bobStrategies=[{
                        y:(choiceOfOutputs[y],comm[y]) for y in range(self.underlyingPolytope.numberOfInputsBob())
                        } 
                            for choiceOfOutputs in 
                                product(*[range(numberOfOutputs) for numberOfOutputs in self.underlyingPolytope.getNumberOfOutputsPerInputBob()])
                            for comm in communicationStrgs]
         
        yield from (Behaviour(self.underlyingPolytope.getBellScenario(),
                        {
                            (inputsAlice,inputsBob,outputsAliceSequence,outputsBob):
                                int(
                                    (outputsBob==stgBob[inputsBob][0])&
                                    (outputsAliceSequence==pairOfStgsAlice[int(stgBob[inputsBob][1])][inputsAlice])
                                    )
                            for (inputsAlice,inputsBob,outputsAliceSequence,outputsBob) in self.underlyingPolytope.getTuplesOfEvents()
                        }) 
                              for stgBob in bobStrategies 
                              for pairOfStgsAlice in aliceStrategiesWithComm
                    )
