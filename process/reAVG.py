#!/usr/bin/env python

from pyAVG.DNAHistoryGraph.graph import DNAHistoryGraph
from pyAVG.inputs.simulator import RandomHistory

def tryExtension(graph):

def reAVG(graph):
	new = copy.copy(graph)
	new.markElem
	while not copy.isAVG():
		copy = tryExtension(copy)
	return copy 

def test_main():
	history = RandomHistory(10, 10)
	avg = history.avg()
	graph = deAVG(avg)
	avg2 = reAVG(graph)
	assert avg2.validate()
	assert avg2.isAVG()
	assert avg2.substitutionCost() >= graph.substitutionCost(lower=True)
	assert avg2.substitutionCost() <= graph.substitutionCost(lower=False)
	assert avg2.rearrangementCost() >= graph.rearrangementCost(lower=True)
	assert avg2.rearrangementCost() <= graph.rearrangementCost(lower=False)

if __name__ == "__main__":
	test_main()
