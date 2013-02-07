#!/usr/bin/env python
import copy
import random

from pyAVG.DNAHistoryGraph.graph import DNAHistoryGraph
from pyAVG.inputs.simulator import RandomHistory

def deAVG(avg, removalDensity=0.0, labelRemovalDensity=0.5, bondRemovalDensity=1.0):
	# Copy graph
	graph = copy.copy(avg)

	# Remove all non-root or leaf segments
	for segment in list(graph.segments):
		if random.random() < removalDensity:
			segment.disconnect()
			graph.segments.remove(segment)
		else:
			if random.random() < labelRemovalDensity:
				segment.label = None
			if len(segment.left.liftedBonds()) > 1:
				segment.left.deleteBond()
			if len(segment.right.liftedBonds()) > 1:
				segment.right.deleteBond()

	# Recompute the event graph from scratch
	graph.eventGraph, graph.segmentThreads = graph.threads()
	graph.timeEventGraph()
	return graph

def test_main():
	for i in range(1000):
		deAVG(RandomHistory(5, 5).avg()).validate()

if __name__ == "__main__":
	test_main()
