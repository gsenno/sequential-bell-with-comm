'''
Created on 9 jul. 2019

@author: gsenno
'''
from bellpolytope import BellPolytope
import cdd as cdd
import numpy as np
import itertools as it
from _functools import reduce
from ncpol2sdpa.sdp_relaxation import imap
import qutip as qt

def getFirstAliceMarginal(nA1,nA2,nB,dist):
    pAlice1 = []
    for x1,y in it.product(range(nA1),range(nB)):
        for a1,b in it.product(range(2),repeat=2):
            pAlice1.append(sum([1/nA2*dist[(b+2*a1+4*a2)+8*(y+nA2*x2+(nA1**2)*x1)]
                for x2 in range(nA2)
                for a2 in range(2)]))
    return pAlice1

def getSecondAliceMarginal(nA1,nA2,nB,dist):
    pAlice2 = []
    for x2,y in it.product(range(nA2),range(nB)):
        for a2,b in it.product(range(2),repeat=2):
            pAlice2.append(sum([1/nA1*dist[(b+2*a1+4*a2)+8*(y+nA2*x2+(nA1**2)*x1)]
                    for x1 in range(nA1)
                    for a1 in range(2)]))
    return pAlice2

def createQubitObservable(unnormalizedBlochVector):
    normalizedBlochVector=1/np.linalg.norm(unnormalizedBlochVector)*np.array(unnormalizedBlochVector)
    paulies=[qt.sigmax(),qt.sigmay(),qt.sigmaz()]
    return sum([paulies[i]*normalizedBlochVector[i] for i in range(0,3)])

def effectForQubitPovm(weight,qubitObservable):
    return weight*(qt.qeye(2)+qubitObservable)

def projectorsForQubitObservable(qubitObservable):
    return list(map(lambda sign : effectForQubitPovm(1/2, sign*qubitObservable),[1,-1]))


def addEffectsForAbortOutcomes(effects,dim):
    return list(map(lambda povm : povm+[0*qt.qeye(dim)],effects))

def createMaxEntState(dimension):
    return 1/np.sqrt(dimension)*(sum([qt.tensor(qt.basis(dimension,i),
                                                   qt.basis(dimension,i)) for i in range(0,dimension)]))

def computeDistributionFromStateAndEffects(state,aliceEffects,bobEffects):
    dist=[]
    for x in range (0,len(aliceEffects)):
        for y in range (0,len(bobEffects)):
            for a in range (0,len(aliceEffects[x])):
                for b in range (0,len(bobEffects[y])):
                    dist.append((qt.tensor(aliceEffects[x][a],bobEffects[y][b])*state*state.dag()).tr())
    return dist                

def createConstraintsMatrixForEfficiencyDual(outputsAliceSequence, outputsBob):
    KAlice = np.sum(outputsAliceSequence)
    KBob = np.sum(outputsBob)
    
    
    totalOutputs = np.sum([a*b for a in outputsAliceSequence for b in outputsBob])
    abortOutputs = totalOutputs - np.sum([(a-1)*(b-1) for a in outputsAliceSequence for b in outputsBob])
    
    NAlice=len(outputsAliceSequence)
    NBob=len(outputsBob)
    totalInputs = NAlice*NBob
     
    nZeroConstraints = 2*abortOutputs
    #nLowerBoundOnDeterministicConstraints = np.prod(outputsAliceSequence)*np.prod(outputsBob)
    nLowerBoundOnDeterministicConstraints = 0
    nUpperBoundOnDeterministicConstraints = np.prod(outputsAliceSequence)*np.prod(outputsBob)
    nConstraints =  nZeroConstraints + nLowerBoundOnDeterministicConstraints + nUpperBoundOnDeterministicConstraints
    nVariables = totalOutputs
 
    matrixHRepresentation=np.zeros((nConstraints,1+nVariables),np.int8)
     
    #IMPOSING THAT THE BELL COEFFICIENTS FOR ABORT EVENTS ARE == 0
    #BOB'S ABORTS
    rowNum=0
    for i in range (0,NAlice):
        for j in range (0,NBob):
            for k in range (0,outputsAliceSequence[i]):
                bellcoeff={}
                bellcoeff[toString(k, outputsBob[j]-1, i, j)]=1
                matrixHRepresentation[rowNum,1:]=toVector(bellcoeff,outputsAliceSequence,outputsBob)
                bellcoeff[toString(k,outputsBob[j]-1,i,j)]=-1
                matrixHRepresentation[rowNum+1,1:]=toVector(bellcoeff,outputsAliceSequence,outputsBob)
                rowNum+=2
            for k in range (0,outputsBob[j]-1):
                bellcoeff={}
                bellcoeff[toString(outputsAliceSequence[i]-1,k,i,j)]=1
                matrixHRepresentation[rowNum,1:]=toVector(bellcoeff,outputsAliceSequence,outputsBob)
                bellcoeff[toString(outputsAliceSequence[i]-1,k,i,j)]=-1
                matrixHRepresentation[rowNum+1,1:]=toVector(bellcoeff,outputsAliceSequence,outputsBob)
                rowNum+=2
    
       
    #IMPOSING THAT THE BELL FUNCTIONAL IS LOWER BOUNDED BY N**2*(N**2-2) ON LOCAL DISTRIBUTIONS WITH ABORT
#     for vector in poly.getVertices():
#         matrixHRepresentation[rowNum,0]=N**2*(N**2-2)
#         matrixHRepresentation[rowNum,1:]=vector
#         rowNum+=1
     
    #IMPOSING THAT THE BELL FUNCTIONAL IS UPPER BOUNDED BY 1 ON LOCAL DISTRIBUTIONS WITH ABORT    
    #vertices=BellPolytope(4,3).getVertices()
    aliceStrategies = generateStrategies(outputsAliceSequence)
    bobStrategies = generateStrategies(outputsBob)
    
    #vertices=[strategyToDistribution(stgAlice,stgBob,outputsAliceSequence,outputsBob) for stgAlice in aliceStrategies for stgBob in bobStrategies]
    for stgAlice in aliceStrategies:
        for stgBob in bobStrategies:
            vertex=strategyToDistribution(stgAlice,stgBob,outputsAliceSequence,outputsBob)
            matrixHRepresentation[rowNum,0]=1
            matrixHRepresentation[rowNum,1:]=-np.array(vertex)
            rowNum+=1
            #matrixHRepresentation[rowNum+1,0]=4**2*(4**2-2)
            #matrixHRepresentation[rowNum+1,1:]=np.array(vector)
    
    return matrixHRepresentation    

def generateLocalVertices(outputsAliceSequence,outputsBob):
    aliceStrategies = generateStrategies(outputsAliceSequence)
    bobStrategies = generateStrategies(outputsBob)
    
    return [strategyToDistribution(stgAlice,stgBob,outputsAliceSequence,outputsBob) for stgAlice in aliceStrategies for stgBob in bobStrategies]
    
    
def toString(a,b,x,y):
    return str(a)+','+str(b)+'|'+str(x)+','+str(y)
    
def toVector(dictionary,outputsAliceSequence,outputsBob):
    vector = []
    for x in range (0,len(outputsAliceSequence)):
        for y in range (0,len(outputsBob)):
            for a in range (0,outputsAliceSequence[x]):
                for b in range (0,outputsBob[y]):
                    if toString(a, b, x, y) in dictionary:
                        vector.append(dictionary[toString(a, b, x, y)])
                    else:
                        vector.append(0)
    return vector
    
def generateStrategies(outputs):
    l=list(it.product(list(range(0,max(outputs))),repeat=len(outputs)))
    return l
    #return list(it.filterfalse(lambda tup : not reduce(lambda x,y : x&y,list(imap(lambda a,b : a<b,tup,outputs))),l))
    
def strategyToDistribution(stgAlice, stgBob,outputsAliceSequence,outputsBob):
    vector = []
    for x in range (0,len(outputsAliceSequence)):
        for y in range (0,len(outputsBob)):
            for a in range (0,outputsAliceSequence[x]):
                for b in range (0,outputsBob[y]):
                    if (a==stgAlice[x])&(b==stgBob[y]):
                        vector.append(1)
                    else:
                        vector.append(0)
#return VerticesToCG(vector, outputsAliceSequence, outputsBob)
    return vector

def VerticesToCG(vector, outputsAliceSequence, outputsBob):
    vertice = []
    #Alice's marginals
    l=0
    for x in range (0,len(outputsAliceSequence)):
        s=0
        for y in range (0,outputsBob[0]):
            s += vector[l+y]
        vertice.append(s)
        for a in range (0,len(outputsBob)):
            l+=outputsAliceSequence[x]*outputsBob[a]
    
    #Bob's marginals
    l=0
    for z in range (0,len(outputsAliceSequence)):
        s=0
        for t in range (0,outputsAliceSequence[0]):
            s += vector[l+outputsBob[z]*t]
        vertice.append(s)
        l+=outputsAliceSequence[0]*outputsBob[z]
        
    #The rest of probabilities
    s=0
    for w in range (0,len(outputsAliceSequence)):
        for x in range (0,len(outputsBob)):
#            vertice.append(vertices[s])
            s+=outputsAliceSequence[w]*outputsBob[x]
    return vertice

def Permutation(vertice,outputsAliceSequence,outputsBob):
    permutedVertice = []
    #Marginals
    for x in range (0, len(outputsBob)):
        permutedVertice.append(vertice[x+len(outputsAliceSequence)])
    for y in range (0, len(outputsAliceSequence)):
        permutedVertice.append(vertice[y])
    
    #The rest of probabilities
    CoefficientMatrix=np.zeros((len(outputsAliceSequence),len(outputsBob)))
    s=0
    #I create a matrix with the rest of probabilities in order to be easy to be permuted
    for l in range (0, len(outputsAliceSequence)):
        for w in range (0, len(outputsBob)):
            CoefficientMatrix[l][w]=vertice[s+len(outputsAliceSequence)+len(outputsBob)]
            s+=1
    
    #I apply the permutation
    for z in range (0, len(outputsAliceSequence)):
        for t in range (0, len(outputsBob)):
            permutedVertice.append(CoefficientMatrix[t][z])
    return permutedVertice

def symmetriseVertices(vertice,permutedVertice):
    symmetricBasis = []
    vector=1/2*(np.array(vertice)+np.array(permutedVertice))
    for element in vector:
        if element not in symmetricBasis:
            symmetricBasis.append(element)
    return symmetricBasis
    
    
def generateVertices1bitOfCommLocalPol(outputsAliceSequence,outputsBob):
    communicationStrgs=list(it.product([0,1],repeat=len(outputsAliceSequence)))
    #strgsAlice = [[[stgAlice[i],comm[i]] for i in range(0,len(stgAlice))] 
    #              for stgAlice in generateStrategies(outputsAliceSequence) for comm in communicationStrgs]
    strgsBob = generateStrategies(list(reduce(lambda acum,elem : acum+[elem,elem],outputsBob,[])))
    
    vertices = []
    for stgAlice in generateStrategies(outputsAliceSequence):
        for comm in communicationStrgs:
            for stgBob in strgsBob:
                vector = []
                for x in range (0,len(outputsAliceSequence)):
                    for y in range (0,len(outputsBob)):
                        for a in range (0,outputsAliceSequence[x]):
                            for b in range (0,outputsBob[y]):
                                if (a==stgAlice[x])&(b==stgBob[2*y+comm[x]]):
                                    vector.append(1)
                                else:
                                    vector.append(0)
                vertices.append(vector)
    return vertices
