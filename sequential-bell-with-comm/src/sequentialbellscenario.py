'''
Created on 16 ago. 2019

@author: gsenno
'''
from bellscenario import BellScenario
from itertools import product

class SequentialBellScenario(BellScenario):


    '''
    numberOfOutputsPerInputAliceSequence is a list of integers such that 
    numberOfOutputsPerInputAlice[0][i]==#outputs for input i of the first Alice and
    equivalently, numberOfOutputsPerInputAlice[1][i] for the second Alice 
    (numberOfOutputsPerInputBob is as in BellScenario).
    '''
    def __init__(self, numberOfOutputsPerInputAliceSequence,numberOfOutputsPerInputBob):
        self.numberOfOutputsPerInputAliceSequence=numberOfOutputsPerInputAliceSequence
        self.outputsAlice1=numberOfOutputsPerInputAliceSequence[0]
        self.outputsAlice2=numberOfOutputsPerInputAliceSequence[1]
        self.numberOfOutputsPerInputBob=numberOfOutputsPerInputBob
        self.tuplesOfEvents = [((x1,x2),y,(a1,a2),b) 
                               for x1 in range(len(self.outputsAlice1))
                               for x2 in range(len(self.outputsAlice2))
                               for y in range(len(self.numberOfOutputsPerInputBob))
                               for a1 in range(self.outputsAlice1[x1])
                               for a2 in range(self.outputsAlice2[x2])
                               for b in range(self.numberOfOutputsPerInputBob[y])]
    
    def numberOfInputsAlice1(self):
        return len(self.outputsAlice1)
    
    def numberOfInputsAlice2(self):
        return len(self.outputsAlice2)
    
    def numberOfInputsAlice(self):
        return self.numberOfInputsAlice1()*self.numberOfInputsAlice2()
    
    def listOfAliceInputs(self):
        return [(x1,x2) for x1,x2 in 
                product(range(self.numberOfInputsAlice1()),range(self.numberOfInputsAlice2()))]
        
    def getNumberOfOutputsPerInputAlice1(self):
        return self.outputsAlice1
    
    def getNumberOfOutputsPerInputAlice2(self):
        return self.outputsAlice2
    
    def isSequentialBehaviour(self,behaviour):
        result=True
        p=behaviour.getProbabilityDictionary()
        for x1 in range(self.numberOfInputsAlice1()):
            for a1 in range(self.getNumberOfOutputsPerInputAlice1()[x1]):
                px1a1 = sum([p[(x1,0),y,(a1,a2),b] for a2 in range(self.getNumberOfOutputsPerInputAlice2()[0])
                            for y in range(self.numberOfInputsBob()) for b in range(self.getNumberOfOutputsPerInputBob()[y])])
                for x2 in range(self.numberOfInputsAlice2()):
                    result = result and px1a1 == sum([p[(x1,x2),y,(a1,a2),b] for a2 in range(self.getNumberOfOutputsPerInputAlice2()[x2])
                            for y in range(self.numberOfInputsBob()) for b in range(self.getNumberOfOutputsPerInputBob()[y])])
        return result
                
                
                
                