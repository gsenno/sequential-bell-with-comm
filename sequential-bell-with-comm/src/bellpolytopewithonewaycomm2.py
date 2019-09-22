'''
Created on 5 ago. 2019

@author: gsenno
'''
from bellpolytope import BellPolytope
from itertools import product
from behaviour import Behaviour

'''
Polytope for bipartite strategies using shared randomness and 1 bit of communication from Alice to Bob.
'''
class BellPolytopeWithOneWayCommunication2(BellPolytope):
    
    def getBellScenario(self):
        return self.underlyingPolytope.getBellScenario()
    
    #see bellpolytope.py for description of params
    def __init__(self, underlyingPolytope):
        self.underlyingPolytope=underlyingPolytope
        
    def getGeneratorForVertices(self):
        #local vertices
        yield from self.underlyingPolytope.getGeneratorForVertices()
        
        communicationStrgs=[format(i,'0'+str(self.underlyingPolytope.numberOfInputsAlice())+'b') 
                            for i in range(1,2**(self.underlyingPolytope.numberOfInputsAlice()-1))]
        listOfAliceInputs=self.underlyingPolytope.listOfAliceInputs() 
        aliceStrategies=[{anInput:(stg[anInput],communication[listOfAliceInputs.index(anInput)]) for anInput in listOfAliceInputs}
                         for stg in self.underlyingPolytope.getAliceStrategies() for communication in communicationStrgs]
        
        bobStrategiesWithComm = [(b1,b2) 
                                   for b1,b2 in product(self.underlyingPolytope.getBobStrategies(),repeat=2)
                                   if not b1==b2]
         
        yield from (Behaviour(self.underlyingPolytope.getBellScenario(),
                        {
                            (inputsAlice,inputsBob,outputsAlice,outputsBob):
                                int(
                                    (outputsAlice==stgAlice[inputsAlice][0])&
                                    (outputsBob==pairOfStgsBob[int(stgAlice[inputsAlice][1])][inputsBob])
                                    )
                            for (inputsAlice,inputsBob,outputsAlice,outputsBob) in self.underlyingPolytope.getTuplesOfEvents()
                        }) 
                              for stgAlice in aliceStrategies 
                              for pairOfStgsBob in bobStrategiesWithComm
                    )
