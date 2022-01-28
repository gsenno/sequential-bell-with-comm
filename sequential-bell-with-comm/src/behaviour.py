'''
Created on 16 ago. 2019

@author: gsenno
'''
from itertools import product

class Behaviour(object):
    
    probabilitiesDictionary={}
    probabilitiesList=[]
    
    def __init__(self, bellScenario, pr):
        self.bellScenario=bellScenario
        if isinstance(pr, list):
            self.probabilitiesList=pr
        else:
            if isinstance(pr, dict):
                self.probabilitiesDictionary=pr
    
    def getProbabilityDictionary(self):
        if self.probabilitiesDictionary=={}:
            i=0
            for event in self.bellScenario.getTuplesOfEvents(): 
                self.probabilitiesDictionary[event]=self.probabilitiesList[i]
                i+=1
        return self.probabilitiesDictionary

    
    def getProbabilityList(self):
        if self.probabilitiesList==[]:
            p=self.getProbabilityDictionary()
            self.probabilitiesList=[p[event] for event in self.bellScenario.getTuplesOfEvents()]
        return self.probabilitiesList
    
    
    def getBipartiteCorrelator(self,inputAlice,inputBob):
        correlator = 0
        probDict = self.getProbabilityDictionary()
        for a,b in product([-1,1],repeat=2):
            correlator = correlator + a*b*probDict[inputAlice,inputBob,a,b]
        return correlator
        
    def getAliceCorrelator(self,inputAlice):
        correlator = 0
        probDict = self.getProbabilityDictionary()
        for a in product([-1,1],repeat=2):
            correlator = correlator + a*b*probDict[inputAlice,inputBob,a,b]
        return correlator
