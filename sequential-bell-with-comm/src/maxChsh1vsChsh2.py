
'''
Created on Feb 12, 2019

@author: gsenno
'''

import numpy as np
from itertools import product

import picos as pic

from bellpolytopewithonewaycomm import BellPolytopeWithOneWayCommunication
from sequentialbellpolytope import SequentialBellPolytope
from sequentialbellscenario import SequentialBellScenario


def CHSH(A,B):
    #chsh bell op give observables
    return pic.kron(A[0],B[0])+pic.kron(A[0],B[1])+pic.kron(A[1],B[0]) \
            -pic.kron(A[1],B[1])

def findQDistMaximizingCHSH2ForCHSH1ValueOf(alpha):
 
    prob=pic.Problem() 
     
    A={}
    for x1,x2 in product(range(2),range(2)):
        for a1,a2 in product(range(2),range(2)):
            A[x1,x2,a1,a2]=prob.add_variable('A_{0}{1}{2}{3}'.format(x1,x2,a1,a2),
                                             (2,2),'hermitian')
            prob.add_constraint(A[x1,x2,a1,a2]>>0)
     
    for x1,x2 in product(range(2),range(2)):
        prob.add_constraint(sum(A[x1,x2,a1,a2] for a1 in range(2) 
                                               for a2 in range(2))==np.eye(2))
    
    #sequential constraints
    for x1,a1 in product(range(2),range(2)):
        prob.add_constraint(sum([A[x1,0,a1,a2] for a2 in (0,1)])==sum([A[x1,1,a1,a2] for a2 in (0,1)]))
    
     
    z0=np.array([[1,0],[0,0]])
    z1=np.array([[0,0],[0,1]])
    x0=np.array([[1,1],[1,1]])/2
    x1=np.array([[1,-1],[-1,1]])/2
     
     
    B={}
    B[0,0]=pic.new_param('B_00',z0)
    B[0,1]=pic.new_param('B_01',z1)
    B[1,0]=pic.new_param('B_10',x0)
    B[1,1]=pic.new_param('B_11',x1)
     
    rho=pic.new_param('rho',np.outer([1,0,0,1],[1,0,0,1])/2)
     
     
    A1=[1/2*sum(A[x1,x2,a1,a2]*(-1)**a1 for a1 in range(2) 
                                        for a2 in range(2) 
                                        for x2 in range(2))
                                                    for x1 in range(2)]
     
    A2=[1/2*sum(A[x1,x2,a1,a2]*(-1)**a2 for a1 in range(2) 
                                        for a2 in range(2) 
                                        for x1 in range(2))
                                                    for x2 in range(2)]
     
    B1=[sum(B[y,b]*(-1)**b for b in range(2)) for y in range(2)]
     
     
    prob.add_constraint(pic.trace(CHSH(A1,B1)*rho)==alpha)
    
    prob.set_objective('max',pic.trace(CHSH(A2,B1)*rho))
     
    prob.solve()
#     
    return [pic.trace(pic.kron(A[x1,x2,a1,a2],B[y,b])*rho).get_value().real
          for x1,x2,y,a1,a2,b in product(range(2),range(2),range(2),range(2),range(2),range(2))] 
    
    

if __name__ == '__main__':
    
    alpha=2.5
    qdist = findQDistMaximizingCHSH2ForCHSH1ValueOf(alpha)
    
     
    outputsAliceSequence = [[2,2],[2,2]]
    outputsBob = [2,2]
    bellScenario = SequentialBellScenario(outputsAliceSequence,outputsBob)
    
    polytope = BellPolytopeWithOneWayCommunication(SequentialBellPolytope(bellScenario))                
    
    if not polytope.contains(qdist):
        print('We found a qdist not reproducible with one bit of comm in the sequential scenario!')