import random
import time
import copy
import os

from pyAVG.DNAHistoryGraph.graph import DNAHistoryGraph
from pyAVG.process.extensionMoves import listCase1, listCase2, listCase3
from pyAVG.utils.tex import *

from pyAVG.inputs.simulator import RandomHistory
from pyAVG.process.deAVG import deAVG
from sonLib.bioio import system
import pyAVG.process.extensionMoves

"""Script generates results for "A Unifying Parsimony Model for Genome Evolution"
in a directory called "results".
"""

def main():
    experimentNumber = 5
    iterationNumber = 1000
    startTime = time.time()
    last = startTime
    results = []
    experiment = 0
    segmentNumber = 5
    epochs = 5
    
    outputDir = "results"
    system("mkdir %s" % outputDir)
    
    while experiment < experimentNumber:
        #Create a random history
        history = RandomHistory(segmentNumber, epochs)
        avg = history.avg()
        
        def breakAllHomologousSides(side, graph):
            if side.bond != None:
                graph.deleteBond(side)
            for child in side.children():
                breakAllHomologousSides(child, graph)
                
        def homologousSidesHaveOnlyTrivialLifts(side):
            if side.bond != None and len(side.nonTrivialLiftedBonds()) > 0:
                return False
            for child in side.children():
                if not homologousSidesHaveOnlyTrivialLifts(child):
                    return False
            return True
        
        sidesToBreak = [ side for side in avg.sides() if side.parent() == None and homologousSidesHaveOnlyTrivialLifts(side) ]
        if len(sidesToBreak) > 0:
            breakAllHomologousSides(sidesToBreak[0], avg)        
        
        #Undo stuff in the first graph
        baseGraph = deAVG(avg)
        assert avg.substitutionAmbiguity() == 0 
        assert avg.rearrangementAmbiguity() == 0
        assert avg.lowerBoundRearrangementCost() == avg.upperBoundRearrangementCost()
        assert baseGraph.lowerBoundRearrangementCost() <= avg.lowerBoundRearrangementCost()
        
        #Selection of histories with what we want
        if avg.lowerBoundRearrangementCost() == 0 or avg.lowerBoundSubstitutionCost() == 0 or baseGraph.substitutionAmbiguity() == 0 or baseGraph.rearrangementAmbiguity() == 0 or len([ segment for segment in baseGraph.segments if len(segment.children) == 0 ]) != 4*segmentNumber:
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
        #assert avg.validate()
        
        #Write stuff about G
        results.append(reportGraph(baseGraph, "G", "n/a", "n/a"))
        #assert baseGraph.validate()
        
        for iteration in range(iterationNumber):
            print "Starting iteration", iteration
            graph = copy.copy(baseGraph)
        
            #Undo the ambiguity
            #lBSC = graph.lowerBoundSubstitutionCost()
            #lBRC = graph.lowerBoundRearrangementCost()
            step = 1
            while graph.ambiguity():
                results.append(reportGraph(graph, "G'", iteration, step))
                c1EL = listCase1(graph)
                c2EL = listCase2(graph)
                c3EL = listCase3(graph)
                #print "There are %s labeling extensions and %s bond extensions" % (len(c1EL), len(c2EL))
                chosenExtension = random.choice(c1EL + c2EL + c3EL)
                chosenExtension.function(chosenExtension.args)
                
                #assert graph.validate()
                #assert lBSC <= graph.lowerBoundSubstitutionCost()
                #assert lBRC <= graph.lowerBoundRearrangementCost()
                #lBSC = graph.lowerBoundSubstitutionCost()
                #lBRC = graph.lowerBoundRearrangementCost()
                
                step += 1
            
            
            #Report final AVG
            for segment in list(graph.segments): #Get rid of useless nodes
                if segment.label == None and segment.left.bond == None and segment.right.bond == None:
                    segment.disconnect()
                    graph.segments.remove(segment) 
        
            # Recompute the event graph from scratch
            graph.eventGraph, graph.segmentThreads = graph.threads()
            graph.timeEventGraph()
            
            
            results.append(reportGraph(graph, "G'", iteration, step))
            #assert graph.validate()
            
            #assert graph.lowerBoundSubstitutionCost() == graph.upperBoundSubstitutionCost()
            #assert graph.lowerBoundRearrangementCost() == graph.upperBoundRearrangementCost()
            #for m in graph.modules():
            #    assert m.isSimple()
                
        experiment += 1
        print 'EXPERIMENT', experiment, time.time() - last
        last = time.time()
    #print "results", results
    
    ###Now format results for printing
    
    #Print table with following fields
    #Experiment, s(H), r(H), lbsc(G), ubsc(G), lbrc(G), ubrc(G), min_G' s(G'),  min_G' r(G'), max_G' s(G'), max_G' r(G'), med_G' s(G'), med_G' s(G')
    
    print "Making summary tables latex"
    
    def fn(statFn, term):
        return int(statFn([ row[term] for row in gPRows ]))
    
    def med(l):
        l = l[:]
        l.sort()
        if len(l) % 2 == 1:
            return l[len(l)/2]
        return (l[len(l)/2] + l[len(l)/2 - 1])/2.0
    
    def getRows():
        experimentResults = [ row for row in results if row["experiment"] == experiment ]
        return [ row for row in experimentResults if row["graphName"] == "H" ][0], [ row for row in experimentResults if row["graphName"] == "G" ][0], [ row for row in experimentResults if row["graphName"] == "G'" and row["ambiguity"] == 0 ]
    
    rTable = [ [ "exp.", "$r(H)$", "$u_r(G)$", "$lbrc(G)$", "$ubrc(G)$", "$r(G_{rmin}')$", "$r(G_{rmax}')$", "$r(G_{rmed}')$" ] ]
    for experiment in range(experimentNumber):
        historyRow, gRow, gPRows = getRows() 
        rTable.append([ str(i) for i in [ experiment, historyRow["lbrc"], gRow["u_r"], gRow["lbrc"], gRow["ubrc"], fn(min, "lbrc"), fn(max, "lbrc"), fn(med, "lbrc") ] ]) 
    
    sTable = [ [ "exp.", "$s(H)$", "$u_s(G)$", "$lbsc(G)$", "$ubsc(G)$", "$s(G_{smin}')$", "$s(G_{smax}')$", "$s(G_{smed}')$" ] ]
    for experiment in range(experimentNumber):
        historyRow, gRow, gPRows = getRows() 
        sTable.append([ str(i) for i in [ experiment, historyRow["lbsc"], gRow["u_s"], gRow["lbsc"], gRow["ubsc"], fn(min, "lbsc"), fn(max, "lbsc"), fn(med, "lbsc")] ]) 

    def writeLatexTable(table, fileName, tableLabel, caption=""):
        fH = open(os.path.join(outputDir, fileName), 'w')
        writeDocumentPreliminaries(fH)
        writePreliminaries(8, fH)
        for line in table:
            writeRow(line, fH)
        writeEnd(fH, tableLabel, caption)
        writeDocumentEnd(fH)
        fH.close()
    
    writeLatexTable(sTable, "aggregateSubs.tex", "subsExpTable", "Results for substitution ambiguity and cost. Starting for an initial evolutionary history H we randomly removed elements to create G and then, by G-bounded extension operations, created a G-bounded AVG G'. Each row represents a separate initial evolutionary history. For each evolutionary history we created 1000 $G$-bounded AVG extensions. $G'_{smin}$, $G'_{smax}$ and $G'_{smax}$ are, respectively, the $G$-bounded extension with minimum, maximum and median substitution cost.")
    writeLatexTable(rTable, "aggregateRearrangements.tex", "rearrangeExpTable", "Follows format of Table \\ref{subsExpTable}.")
        
    def writeTSVTable(table, fileName):
        fH = open(os.path.join(outputDir, fileName), 'w')
        for line in table:
            fH.write("\t".join(line) + "\n")
        fH.close()
        
    writeTSVTable(sTable, "aggregateSubs.tsv")
    writeTSVTable(rTable, "aggregateRearrangements.tsv")
    
    print "Making step-wise .tsv files"
    
    #Write .tsv files showing the change in the bounds during extension from G to G' during an iteration
    
    def writeStepFile(argName, fn):
        fH = open(os.path.join(outputDir, "%s.%s.steps.tsv" % (experiment, argName)), 'w')
        historyRow = getRows()[0]
        fH.write("#TSV showing %s value, original history had %s value \n" % (argName, fn(historyRow)))
        experimentResults = [ row for row in results if row["experiment"] == experiment and row["graphName"] == "G'" ]
        experimentResultsByIteration =  [ [] for i in range(iterationNumber) ]
        for row in experimentResults:
            experimentResultsByIteration[row["iteration"]].append(row)
        #Print row for scale
        fH.write("\t".join([ str(j) for j in range(max([ len(i) for i in experimentResultsByIteration ])) ]) + "\n")
        for it in xrange(iterationNumber):
            itRows = experimentResultsByIteration[it]
            fH.write("%s\n" % "\t".join([ str(fn(row)) for row in itRows]))
        fH.close()
    
    for experiment in range(experimentNumber):
        writeStepFile("lbrc", lambda row : row["lbrc"])
        writeStepFile("ubrc", lambda row : row["ubrc"])
        writeStepFile("ubrc-lbrc", lambda row : row["ubrc"] - row["lbrc"])
        writeStepFile("lbsc", lambda row : row["lbsc"])
        writeStepFile("ubsc", lambda row : row["ubsc"])
        writeStepFile("lbsc-ubsc", lambda row : row["ubsc"] - row["lbsc"])

    print "Making graphviz plots"

    #Write the dot files of the H, G and G_min, G_max and G_med for each experiment
    for experiment in range(experimentNumber):
        historyRow, gRow, gPRows = getRows() 
        dirName = os.path.join(outputDir, "%s_graphViz" % experiment)
        system("mkdir %s" % dirName)
        def fn(statFn, term):
            i = statFn([ row[term] for row in gPRows ])
            return [ row for row in gPRows if int(row[term]) == i ][0]
        def writeDot(fileName, row):
            fileName = os.path.join(dirName, fileName)
            fH = open(fileName, 'w')
            fH.write("%s\n" % row["dot"])
            fH.close()
            system("dot %s -Tpdf > %s.pdf" % (fileName, fileName))
        writeDot("history", historyRow)
        writeDot("g", gRow)
        writeDot("gPRRMin" , fn(min, "ubrc"))
        writeDot("gPRRMax", fn(max, "ubrc"))
        writeDot("gPRSMin" , fn(min, "ubsc"))
        writeDot("gPRSMax", fn(max, "ubsc"))
        
    print "Simulations took %s seconds" % (time.time() - startTime)

if __name__ == '__main__':
    main()