#!/usr/bin/env python

from pyAVG.inputs.simulator import RandomHistory
from deAVG import deAVG
import moves

def listPossibleExtensions(graph):
	return sum([detector(graph) for detector in moves.detectors], [])

def tryExtension(graph):
	possibleExtensions = listPossibleExtensions(graph)
	assert len(possibleExtensions) > 0
	chosenExtension = random.random(possibleExtensions)	
	return chosenExtension.function(chosenExtension.data)

def graphElements(graph):
	return list(graph.segments) + filter(lambda X: X.bond is not None, graph.sides()) + [segment.label for segment in graph.segments if segment.label is not None]

def reAVG(graph):
	new = copy.copy(graph)
	new.irreducibleElements = graphElements(new)
	while not new.isAVG():
		new = tryExtension(new)
	return new 

def test_main():
	history = RandomHistory(10, 10)
	graph = deAVG(history.avg())
	avg2 = reAVG(graph)
	assert avg2.validate()
	assert avg2.isAVG()
	assert avg2.substitutionCost() >= graph.substitutionCost(lower=True)
	assert avg2.substitutionCost() <= graph.substitutionCost(lower=False)
	assert avg2.rearrangementCost() >= graph.rearrangementCost(lower=True)
	assert avg2.rearrangementCost() <= graph.rearrangementCost(lower=False)

if __name__ == "__main__":
	test_main()
