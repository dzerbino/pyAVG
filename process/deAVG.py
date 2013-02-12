#!/usr/bin/env python
import copy
import random

from pyAVG.DNAHistoryGraph.graph import DNAHistoryGraph
from pyAVG.inputs.simulator import RandomHistory

def deAVG(avg, removalDensity=0.1, labelRemovalDensity=0.5, bondRemovalDensity=0.25):
	# Copy graph
	graph = copy.copy(avg)

	# Remove all non-root or leaf segments
	for segment in list(graph.segments):
		if random.random() < removalDensity and segment.parent != None:
			segment.disconnect()
			graph.segments.remove(segment)
		else:
			if random.random() < labelRemovalDensity:
				segment.label = None
			if random.random() < bondRemovalDensity:
				segment.left.deleteBond()
			if random.random() < bondRemovalDensity:
				segment.right.deleteBond()
	
	for segment in list(graph.segments): #Get rid of useless nodes
		if segment.parent != None and segment.label == None and segment.left.bond == None and segment.right.bond == None and len(segment.children) <= 1:
			segment.disconnect()
			graph.segments.remove(segment) 

	# Recompute the event graph from scratch
	graph.eventGraph, graph.segmentThreads = graph.threads()
	graph.timeEventGraph()
	return graph

def test_main():
	for i in range(1000):
		deAVG(RandomHistory(5, 5).avg()).validate()

if __name__ == "__main__":
	test_main()
