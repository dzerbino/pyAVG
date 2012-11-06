#!/usr/bin/env python

import random
import time

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
	assert new.validate()
	count = 0
	while not new.isAVG():
		tryExtension(new)
		assert new.validate()
		new.makeGBounded()
		assert new.validate()
		count += 1
		if count > 1000:
			# Poor man's infinite loop trap
			print new.dot()
			assert False
	return new 

def test_main():
	last = time.time()
	for i in range(1000):
		print 'EXPERIMENT', i, time.time() - last
		last = time.time()
		history = RandomHistory(5, 5)
		avg = history.avg()
		assert avg.validate()
		graph = deAVG(avg)
		assert graph.validate()
		avg2 = reAVG(graph)
		assert avg2.validate()
		assert avg2.isAVG()
		assert avg2.substitutionCost() >= graph.substitutionCost()
		assert avg2.rearrangementCost() >= graph.rearrangementCost()

if __name__ == "__main__":
	test_main()
