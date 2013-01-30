import unittest
import random
import time

from pyAVG.DNAHistoryGraph.graph import DNAHistoryGraph
from pyAVG.process.extensionMoves import listCase1

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
        
    def testCase1_random(self):
        last = time.time()
        for i in range(100):
            print 'EXPERIMENT', i, time.time() - last
            last = time.time()
            history = RandomHistory(10, 10)
            avg = history.avg()
            assert avg.validate()
            graph = deAVG(avg)
            assert graph.validate()
            print "Avg has substitution ambiguity %s, lbsc %i and ubsc %i" % (avg.substitutionAmbiguity(), avg.lowerBoundSubstitutionCost(), avg.upperBoundSubstitutionCost())
            while graph.substitutionAmbiguity():
                print "Graph has substitution ambiguity %s, lbsc %i and ubsc %i" % (graph.substitutionAmbiguity(), graph.lowerBoundSubstitutionCost(), graph.upperBoundSubstitutionCost()) 
                chosenExtension = random.choice(listCase1(graph))
                chosenExtension.function(chosenExtension.args)
            print "Finally graph has substitution ambiguity %s, lbsc %i and ubsc %i" % (graph.substitutionAmbiguity(), graph.lowerBoundSubstitutionCost(), graph.upperBoundSubstitutionCost()) 
            print graph.dot()
if __name__ == '__main__':
    unittest.main()