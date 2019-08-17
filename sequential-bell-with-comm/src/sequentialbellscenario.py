'''
Created on 16 ago. 2019

@author: gsenno
'''
from bellscenario import BellScenario

class SequentialBellScenario(BellScenario):
    '''
    classdocs
    '''


    def __init__(self, outputsAliceSequence,numberOfOutputsPerInputBob):
        self.outputsAliceSequence=outputsAliceSequence
        self.outputsAlice1=outputsAliceSequence[0]
        self.outputsAlice2=outputsAliceSequence[1]
        self.numberOfOutputsPerInputBob=numberOfOutputsPerInputBob
        self.tuplesOfEvents = [((x1,x2),y,(a1,a2),b) 
                               for x1 in range(len(self.outputsAlice1))
                               for x2 in range(len(self.outputsAlice2))
                               for y in range(len(self.numberOfOutputsPerInputBob))
                               for a1 in range(self.outputsAlice1[x1])
                               for a2 in range(self.outputsAlice1[x2])
                               for b in range(self.numberOfOutputsPerInputBob[y])]
    
    def numberOfInputsAlice1(self):
        return len(self.outputsAlice1)
    
    def numberOfInputsAlice2(self):
        return len(self.outputsAlice2)
    
    def getNumberOfOutputsPerInputAlice1(self):
        return self.outputsAlice1
    
    def getNumberOfOutputsPerInputAlice2(self):
        return self.outputsAlice2