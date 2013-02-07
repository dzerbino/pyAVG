import unittest
import random
import time

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
            history = RandomHistory(10, 10)
            avg = history.avg()
            def writeGraph(graph, file):
                fileHandle = open(file, 'w')
                fileHandle.write("%s\n" % avg.dot())
                fileHandle.close()
                system("dot -Tpdf %s > %s.pdf" % (file, file))
            writeGraph(avg, "history.dot")
            print "Avg has substitution ambiguity %s, lbsc %i and ubsc %i" % (avg.substitutionAmbiguity(), avg.lowerBoundSubstitutionCost(), avg.upperBoundSubstitutionCost())
            assert avg.validate()
            graph = deAVG(avg)
            assert graph.validate()
            writeGraph(graph, "graph.dot")
            i = graph.substitutionAmbiguity()
            while graph.substitutionAmbiguity():
                print "Graph has substitution ambiguity %s, lbsc %i and ubsc %i" % (graph.substitutionAmbiguity(), graph.lowerBoundSubstitutionCost(), graph.upperBoundSubstitutionCost()) 
                chosenExtension = random.choice(listCase1(graph))
                chosenExtension.function(chosenExtension.args)
            
            while graph.rearrangementAmbiguity():
                print "Graph has rearrangement ambiguity %s, lbsc %i and ubsc %i" % (graph.rearrangementAmbiguity(), graph.lowerBoundSubstitutionCost(), graph.upperBoundSubstitutionCost()) 
                chosenExtension = random.choice(listCase2(graph))
                chosenExtension.function(chosenExtension.args)
                
            print "Finally graph has substitution ambiguity %s, lbsc %i and ubsc %i" % (graph.substitutionAmbiguity(), graph.lowerBoundSubstitutionCost(), graph.upperBoundSubstitutionCost()) 
            print "hello", i
            writeGraph(graph, "avg.dot")

if __name__ == '__main__':
    unittest.main()