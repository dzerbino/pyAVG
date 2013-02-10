import unittest
import random
import time
import copy
import sys

from pyAVG.DNAHistoryGraph.graph import DNAHistoryGraph
from pyAVG.process.extensionMoves import listCase1, listCase2

from pyAVG.inputs.simulator import RandomHistory
from deAVG import deAVG
from sonLib.bioio import system
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
        
    def testCase1_random(self):
        last = time.time()
        for i in range(1):
            print 'EXPERIMENT', i, time.time() - last
            last = time.time()
            
            #Create a random history
            history = RandomHistory(3, 3)
	    print history.dot()
            avg = history.avg()
            
            #Functions for reporting the results
            def writeGraph(graph, file):
                graph = copy.copy(graph)
                graph.addFreeRoots()
                fileHandle = open(file, 'w')
                fileHandle.write("%s\n" % graph.dot())
                fileHandle.close()
                system("dot -Tjpg %s > %s.jpg" % (file, file))
            
            def reportGraph(graph, graphName):
                print "%s has u %s, u_s %s, u_r %s, lbsc %i, ubsc %i, lbrc %i, ubrc %i" % (graphName, \
                                                                                           graph.ambiguity(), \
                                                                                           graph.substitutionAmbiguity(), \
                                                                                           graph.rearrangementAmbiguity(), \
                                                                                           graph.lowerBoundSubstitutionCost(), \
                                                                                           graph.upperBoundSubstitutionCost(), \
                                                                                           graph.lowerBoundRearrangementCost(), \
                                                                                           graph.upperBoundRearrangementCost())
            
            #Report the starting point
            reportGraph(avg, "AVG")
            writeGraph(avg, "history.dot")
            assert avg.validate()
	    sys.exit(0)
            
            #Undo stuff in the first graph
            graph = deAVG(avg)
            
            #Write stuff about G
            writeGraph(graph, "graph.dot")
            reportGraph(graph, "G")
            assert graph.validate()
            
            #Undo the ambiguity
            lBSC = graph.lowerBoundSubstitutionCost()
            lBRC = graph.lowerBoundRearrangementCost()
            while graph.ambiguity():
                c1EL = listCase1(graph)
                c2EL = listCase2(graph)
                print "There are %s labeling extensions and %s bond extensions" % (len(c1EL), len(c2EL))
                chosenExtension = random.choice(c1EL + c2EL)
                chosenExtension.function(chosenExtension.args)
                
                reportGraph(graph, "G'")
                assert graph.validate()
                assert lBSC <= graph.lowerBoundSubstitutionCost()
                assert lBRC <= graph.lowerBoundRearrangementCost()
                lBSC = graph.lowerBoundSubstitutionCost()
                lBRC = graph.lowerBoundRearrangementCost()
            
            #Report final AVG
            reportGraph(graph, "H")
            writeGraph(graph, "avg.dot")
            assert graph.validate()
            
            assert graph.lowerBoundSubstitutionCost() == graph.upperBoundSubstitutionCost()
            assert graph.lowerBoundRearrangementCost() == graph.upperBoundRearrangementCost()
            
            for m in graph.modules():
                assert m.isSimple()

if __name__ == '__main__':
    unittest.main()
