import unittest
import random
import time
import copy

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
        
    def testCases_random(self):
        experimentNumber = 10
        iterationNumber = 10
        last = time.time()
        results = []
        experiment = 0
        while experiment < experimentNumber:
            print 'EXPERIMENT', experiment, time.time() - last
            last = time.time()
            
            #Create a random history
            history = RandomHistory(10, 3)
            avg = history.avg()
            
            #Undo stuff in the first graph
            baseGraph = deAVG(avg)
            
            if baseGraph.substitutionAmbiguity() == 0 or baseGraph.rearrangementAmbiguity() == 0:
                continue
            
            def reportGraph(graph, graphName, iteration, step):
                graph = copy.copy(graph)
                graph.addFreeRoots()
                return { "graphName":graphName, 
                                "experiment":experiment,
                                "iteration":iteration, 
                                "step":step, "ambiguity":graph.ambiguity(),
                                "u_s":graph.substitutionAmbiguity(),
                                "u_r":graph.rearrangementAmbiguity(),
                                 "lbsc":graph.lowerBoundSubstitutionCost(),
                                 "lbrc":int(graph.lowerBoundRearrangementCost()),
                                 "ubsc":graph.upperBoundSubstitutionCost(),
                                 "ubrc":int(graph.upperBoundRearrangementCost()),
                                 "dot":("%s\n" % graph.dot()) }
            
            #Report the starting point
            results.append(reportGraph(avg, "H", "n/a", "n/a"))
            assert avg.validate()
            
            #Write stuff about G
            results.append(reportGraph(baseGraph, "G", "n/a", "n/a"))
            assert baseGraph.validate()
            
            for iteration in range(iterationNumber):
                print "Starting iteration", iteration
                graph = copy.copy(baseGraph)
            
                #Undo the ambiguity
                lBSC = graph.lowerBoundSubstitutionCost()
                lBRC = graph.lowerBoundRearrangementCost()
                step = 1
                while graph.ambiguity():
                    results.append(reportGraph(graph, "G'", iteration, step))
                    c1EL = listCase1(graph)
                    c2EL = listCase2(graph)
                    #print "There are %s labeling extensions and %s bond extensions" % (len(c1EL), len(c2EL))
                    chosenExtension = random.choice(c1EL + c2EL)
                    chosenExtension.function(chosenExtension.args)
                    
                    #assert graph.validate()
                    #assert lBSC <= graph.lowerBoundSubstitutionCost()
                    #assert lBRC <= graph.lowerBoundRearrangementCost()
                    #lBSC = graph.lowerBoundSubstitutionCost()
                    #lBRC = graph.lowerBoundRearrangementCost()
                    
                    step += 1
                
                #Report final AVG
                results.append(reportGraph(graph, "G'", iteration, step))
                assert graph.validate()
                
                assert graph.lowerBoundSubstitutionCost() == graph.upperBoundSubstitutionCost()
                assert graph.lowerBoundRearrangementCost() == graph.upperBoundRearrangementCost()
                for m in graph.modules():
                    assert m.isSimple()
                    
            experiment += 1
        #print "results", results
        
        ###Now format results for printing
        
        #Print table with following fields
        #Experiment, s(H), r(H), lbsc(G), ubsc(G), lbrc(G), ubrc(G), min_G' s(G'),  min_G' r(G'), max_G' s(G'), max_G' r(G'), med_G' s(G'), med_G' s(G')
        
        print "Summary tables latex"
        
        def fn(statFn, term):
            return statFn([ row[term] for row in gPRows ])
        
        def med(l):
            l = l[:]
            l.sort()
            if len(l) % 2 == 1:
                return l[len(l)/2]
            return (l[len(l)/2] + l[len(l)/2 - 1])/2.0
        
        def getRows():
            experimentResults = [ row for row in results if row["experiment"] == experiment ]
            return [ row for row in experimentResults if row["graphName"] == "H" ][0], [ row for row in experimentResults if row["graphName"] == "G" ][0], [ row for row in experimentResults if row["graphName"] == "G'" and row["ambiguity"] == 0 ]
        
        fH = open("aggregateSubs.tex", 'w')
        print "\t".join([ "experiment", "r(H)", "u_r(G)", "lbrc(G)", "ubrc(G)", "min_G' r(G')", "max_{G'} r(G')", "med_G'r(G')" ])
        for experiment in range(experimentNumber):
            historyRow, gRow, gPRows = getRows() 
            print "\t".join([ str(i) for i in [ experiment, historyRow["lbrc"], gRow["u_r"], gRow["lbrc"], gRow["ubrc"], fn(min, "lbrc"), fn(max, "lbrc"), fn(med, "lbrc") ] ]) 
        
        print "\t".join([ "experiment", "s(H)", "u_s(G)", "lbsc(G)", "ubsc(G)", "s_{min}(G')", "max_G' s(G')", "med_G' s(G')", ])
        for experiment in range(experimentNumber):
            historyRow, gRow, gPRows = getRows() 
            print "\t".join([ str(i) for i in [ experiment, historyRow["lbsc"], gRow["u_s"], gRow["lbsc"], gRow["ubsc"], fn(min, "lbsc"), fn(max, "lbsc"), fn(med, "lbsc")] ]) 

        

if __name__ == '__main__':
    unittest.main()