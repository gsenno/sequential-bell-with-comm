import numpy as np
from bellpolytope import BellPolytope
from bellscenario import BellScenario
from bellpolytopewithonewaycomm import BellPolytopeWithOneWayCommunication
from ppl import Variable, Generator_System, C_Polyhedron, point
from sequentialbellpolytope import SequentialBellPolytope
from sequentialbellscenario import SequentialBellScenario

if __name__ == '__main__':

    outputsAliceSeq = [[2,2],[2,2]]
    outputsBob = [2,2]
    scenario = SequentialBellScenario(outputsAliceSeq,outputsBob)
    
    variables=[Variable(i) for i in range(len(scenario.getTuplesOfEvents()))]
    gs = Generator_System()
    
    for v in BellPolytopeWithOneWayCommunication(SequentialBellPolytope(scenario)).getGeneratorForVertices():
        prob=v.getProbabilityList()
        gs.insert(point(sum(prob[i]*variables[i] for i in range(len(prob)))))
        
    poly = C_Polyhedron(gs)
    constraints = poly.constraints()
    
    for constraint in constraints:
        inequality = str(constraint.inhomogeneous_term().__float__())+' '
        for coef in constraint.coefficients():
            inequality = inequality+str(-coef.__float__())+' '
        print(inequality)
    
    
