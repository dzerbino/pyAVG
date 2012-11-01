#!/usr/bin/env python
import copy
import random

from pyAVG.DNAHistoryGraph.graph import DNAHistoryGraph
from pyAVG.inputs.simulator import RandomHistory

def deAVG(avg, removalDensity=0.5):
	# Copy graph
	graph = copy.copy(avg)

	# Remove all non-root or leaf segments
	for segment in list(graph.segments):
		if random.random() < removalDensity:
			segment.disconnect()
			graph.segments.remove(segment)

	# Recompute the event graph from scratch
	graph.eventGraph, graph.segmentThreads = graph.threads()
	graph.timeEventGraph()
	return graph

def test_main():
	history = RandomHistory(5, 5)
	avg = history.avg()
	graph = deAVG(avg)
	assert graph.validate()

if __name__ == "__main__":
	test_main()
