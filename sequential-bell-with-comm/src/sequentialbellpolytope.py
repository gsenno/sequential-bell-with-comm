'''
Created on 8 ago. 2019

@author: gsenno
'''
from itertools import product
from _functools import reduce
from bellpolytope import BellPolytope

class SequentialBellPolytope(BellPolytope):
    '''
    classdocs
    '''
    def __init__(self, sequentialBellScenario):
        BellPolytope.__init__(self, sequentialBellScenario)
        
    def getAliceStrategies(self):
        alice1strategies = [{x1:choiceOfOutputs[x1] for x1 in range(self.bellScenario.numberOfInputsAlice1())}
                            for choiceOfOutputs in 
                                product(*[range(numberOfOutputs) for numberOfOutputs in self.bellScenario.getNumberOfOutputsPerInputAlice1()])]
        
        extendOutputsAlice2 = reduce(lambda x,y : x+y,[self.bellScenario.numberOfInputsAlice1()*[x] for x in self.bellScenario.getNumberOfOutputsPerInputAlice2()],[])
        alice2strategies = [
                                {
                                (x1,x2):choiceOfOutputs[self.bellScenario.numberOfInputsAlice1()*x1+x2] 
                                for x1,x2 in product(range(self.bellScenario.numberOfInputsAlice1()),range(self.bellScenario.numberOfInputsAlice2()))
                                }
                            for choiceOfOutputs in product(*[range(numberOfOutputs) for numberOfOutputs in extendOutputsAlice2])
                            ]
        
        return [
                {
                (x1,x2):(stgAlice1[x1],stgAlice2[(x1,x2)])
                for x1,x2 in product(range(self.bellScenario.numberOfInputsAlice1()),range(self.bellScenario.numberOfInputsAlice2()))
                }
                for stgAlice1 in alice1strategies for stgAlice2 in alice2strategies]