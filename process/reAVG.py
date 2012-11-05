#!/usr/bin/env python

import random

from pyAVG.inputs.simulator import RandomHistory
from pyAVG.DNAHistoryGraph.extension import GraphExtension
from pyAVG.DNAHistoryGraph.graph import DNAHistoryGraph
from deAVG import deAVG
import extensionMoves

def listPossibleExtensions(graph):
	return sum([detector(graph) for detector in extensionMoves.detectors], [])

def tryExtension(graph):
	possibleExtensions = listPossibleExtensions(graph)
	assert len(possibleExtensions) > 0
	chosenExtension = random.choice(possibleExtensions)	
	return chosenExtension.function(chosenExtension.args)

def reAVG(graph):
	new = GraphExtension(graph)
	count = 0
	while not new.isAVG():
		tryExtension(new)
		new.makeGBounded()
		count += 1
		if count > 1000:
			print new.dot()
			assert False
	return new 

def test_main():
	history = RandomHistory(5, 5)
	graph = deAVG(history.avg())
	avg2 = reAVG(graph)
	assert avg2.validate()
	assert avg2.isAVG()
	assert avg2.substitutionCost() >= graph.substitutionCost(lowerBound=True)
	assert avg2.rearrangementCost() >= graph.rearrangementCost(lowerBound=True)

if __name__ == "__main__":
	test_main()
