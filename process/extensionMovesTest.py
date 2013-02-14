import unittest
import random
import time
import copy
import os

from pyAVG.DNAHistoryGraph.graph import DNAHistoryGraph
from pyAVG.process.extensionMoves import listCase1, listCase2
from pyAVG.utils.tex import *

from pyAVG.inputs.simulator import RandomHistory
from deAVG import deAVG
import extensionMoves

class ExtensionMovesTest(unittest.TestCase):
    
    def setUp(self):
        unittest.TestCase.setUp(self)
        self.g = DNAHistoryGraph()
        self.s1 = self.g.newSegment()
        self.s2 = self.g.newSegment("A")
        self.s3 = self.g.newSegment("T")
        self.s4 = self.g.newSegment("A")
        self.g.createBranch(self.s1, self.s2)
        self.g.createBranch(self.s1, self.s3)
        self.g.createBranch(self.s1, self.s4)
    
    def tearDown(self):
        unittest.TestCase.tearDown(self)
        
    def testCase1(self):
        print "Graph has substitution ambiguity", self.g.substitutionAmbiguity()
        while self.g.substitutionAmbiguity():
            chosenExtension = random.choice(listCase1(self.g))
            chosenExtension.function(chosenExtension.args)
            print "Graph now has substitution ambiguity", self.g.substitutionAmbiguity()
        self.assertEquals(self.g.substitutionAmbiguity(), 0)
        
    def testCases_random(self):
        experimentNumber = 10
        iterationNumber = 10
        last = time.time()
        experiment = 0
        while experiment < experimentNumber:
            print 'EXPERIMENT', experiment, time.time() - last
            last = time.time()
            
            #Create a random history
            history = RandomHistory(3, 3)
            avg = history.avg()
            
            #Undo stuff in the first graph
            baseGraph = deAVG(avg)
            
            if baseGraph.substitutionAmbiguity() == 0 or baseGraph.rearrangementAmbiguity() == 0:
                continue
            
            def reportGraph(graph, graphName, iteration, step):
                graph = copy.copy(graph)
                graph.addFreeRoots()
                print "\t".join([ "graphName", graphName,
                                "experiment", str(experiment),
                                "iteration", str(iteration), 
                                "step", str(step), 
                                "ambiguity", str(graph.ambiguity()),
                                "u_s", str(graph.substitutionAmbiguity()),
                                "u_r", str(graph.rearrangementAmbiguity()),
                                 "lbsc", str(graph.lowerBoundSubstitutionCost()),
                                 "lbrc", str(graph.lowerBoundRearrangementCost()),
                                 "ubsc", str(graph.upperBoundSubstitutionCost()),
                                 "ubrc", str(graph.upperBoundRearrangementCost()) ])
            
            #Report the starting point
            reportGraph(avg, "H", "n/a", "n/a")
            assert avg.validate()
            
            #Write stuff about G
            reportGraph(baseGraph, "G", "n/a", "n/a")
            assert baseGraph.validate()
            
            for iteration in range(iterationNumber):
                print "Starting iteration", iteration
                graph = copy.copy(baseGraph)
            
                #Undo the ambiguity
                lBSC = graph.lowerBoundSubstitutionCost()
                lBRC = graph.lowerBoundRearrangementCost()
                step = 1
                while graph.ambiguity():
                    reportGraph(graph, "G'", iteration, step)
                    c1EL = listCase1(graph)
                    c2EL = listCase2(graph)
                    #print "There are %s labeling extensions and %s bond extensions" % (len(c1EL), len(c2EL))
                    chosenExtension = random.choice(c1EL + c2EL)
                    chosenExtension.function(chosenExtension.args)
                    
                    assert graph.validate()
                    assert lBSC <= graph.lowerBoundSubstitutionCost()
                    assert lBRC <= graph.lowerBoundRearrangementCost()
                    lBSC = graph.lowerBoundSubstitutionCost()
                    lBRC = graph.lowerBoundRearrangementCost()
                    
                    step += 1
                
                #Report final AVG
                reportGraph(graph, "G'", iteration, step)
                assert graph.validate()
                
                assert graph.lowerBoundSubstitutionCost() == graph.upperBoundSubstitutionCost()
                assert graph.lowerBoundRearrangementCost() == graph.upperBoundRearrangementCost()
                for m in graph.modules():
                    assert m.isSimple()
                    
            experiment += 1
       

if __name__ == '__main__':
    unittest.main()