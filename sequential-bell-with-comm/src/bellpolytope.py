'''
Created on 27 dic. 2018

@author: gsenno
'''
from mosek.fusion import *
import numpy as np
import picos as pic
from itertools import product
from behaviour import Behaviour

class BellPolytope:
    
    # numberOfOutputsPerInputAlice is a list of integers such that numberOfOutputsPerInputAlice[i]==#outputs for Alice's input i (idem numberOfOutputsPerInputBob).
    # the lists do not need to be of the same length.
    def __init__(self,bellScenario):
        self.bellScenario = bellScenario
    
    def getBellScenario(self):
        return self.bellScenario
    
    def getTuplesOfEvents(self):
        return self.bellScenario.getTuplesOfEvents()
    
    def getNumberOfOutputsPerInputAlice(self):
        return self.bellScenario.getNumberOfOutputsPerInputAlice()
    
    def getNumberOfOutputsPerInputBob(self):
        return self.bellScenario.getNumberOfOutputsPerInputBob()
    
    def numberOfInputsAlice(self):
        return self.bellScenario.numberOfInputsAlice()
    
    def numberOfInputsBob(self):
        return self.bellScenario.numberOfInputsBob()
    
    def _strategiesGenerator(self,outputs):
        yield from ({anInput:choiceOfOutputs[anInput] for anInput in range(len(outputs))}
                            for choiceOfOutputs in product(*[range(numberOfOutputs) for numberOfOutputs in outputs]))
        
    def getAliceStrategies(self):
        return self._strategiesGenerator(self.getNumberOfOutputsPerInputAlice())
    
    def getBobStrategies(self):
        return self._strategiesGenerator(self.getNumberOfOutputsPerInputBob())

    def _strategyToBehaviour(self,stgAlice, stgBob):
        distribution = {(inputsAlice,inputsBob,outputsAliceSequence,outputsBob):
                            int(
                                (outputsAliceSequence==stgAlice[inputsAlice])&
                                (outputsBob==stgBob[inputsBob])
                                )
                        for (inputsAlice,inputsBob,outputsAliceSequence,outputsBob) in self.getTuplesOfEvents()}
        return Behaviour(self.bellScenario,distribution)
    
    def getGeneratorForVertices(self):
        return (self._strategyToBehaviour(stgAlice, stgBob)
                 for stgAlice in self.getAliceStrategies()
                 for stgBob in self.getBobStrategies()) 
    
    def getListOfVertices(self):
        return list(self.getGeneratorForVertices())

    def contains(self,distribution):
    # Tests if the point distribution is inside the convex Hull of the points D
    # distribution should be a np multidim array
    # D should be any np multidim array with first index labelling the points of convex set
    # output is list containing the solver status and the solution 
    #reshape so we have vectors
        D=list(map(lambda vertice : vertice.getProbabilityList(), self.getListOfVertices()))
        N=len(D)
        D=pic.new_param('D',D)
        #define problem
        prob=pic.Problem()
        #cerate prob vector
        p=prob.add_variable('p',N)
        #add desired point
        distribution=np.reshape(distribution,[1,-1])
        distribution=pic.new_param('distribution',distribution)
        #feasibilitiy test
        prob.set_objective('min',0*p[0])
         
        #constraints: positivity, normalisation, correct vector
        prob.add_constraint(p>=0)
        prob.add_constraint([[1 for __ in range(N)]]*p==1)
    #    prob.add_constraint(pic.sum([p[i] for i in range(N)])==1.0)
        prob.add_constraint(p.T*D==distribution)
     
        prob.solve(verbose=1)
    #    print(prob)
         
        #get optimal variables and reshape
        pvalues=np.array(p.value)
        return prob.status=='optimal'
        #return [prob.status,pvalues]

#     def contains(self,dist):
#         with Model("lo1") as M:
#             
#             vertices=list(map(lambda vertice : vertice.getProbabilityList(), self.getListOfVertices()))
#             # Create variables
#             bellFunctional = M.variable("func",len(vertices[0]))
#             localBound = M.variable("bound", 1)
#     
#             # Create constraints
#             vertexNum=0
#             for vertex in vertices:
#                 M.constraint('const'+str(vertexNum),Expr.sub(Expr.dot(bellFunctional,vertex),localBound)
#                              ,Domain.lessThan(0))
#                 vertexNum+=1
#                 
#             M.constraint('norm',Expr.sub(Expr.dot(bellFunctional,dist),localBound),Domain.lessThan(1))
#             
#             # Set the objective function to (c^t * x)
#             M.objective("obj", ObjectiveSense.Maximize, Expr.sub(Expr.dot(bellFunctional,dist),localBound))
#     
#             # Solve the problem
#             M.solve()
#             result = [M.primalObjValue(),bellFunctional.level(),localBound.level()]
#             
#             return result[0]<=result[2]
