#!/usr/bin/env python

class Extension(object):
	def __init__(self, data, function):
		self.args = args
		self.function = function

def select_oldest_segments(selected, segment, graph):
	if segment in graph.irreducible:
		return selected

	for other in list(selected):
		if not graph.eventGraph.testConstraint(graph.segmentThreads[segment], graph.segmentThreads[other]):
			return selected
		elif not graph.eventGraph.testConstraint(graph.segmentThreads[other], graph.segmentThreads[segment]):
			selected.remove(other)
	return selected

def oldest_segments(graph):
	return reduce(lambda X, Y: select_oldest_segments(X, Y, graph), graph.segments, [])

def select_oldest_sides(selected, side, graph):
	if side in graph.irreducible:
		return selected

	for other in list(selected):
		if not graph.eventGraph.testConstraint(graph.segmentThreads[side.segment], graph.segmentThreads[other.segment]):
			return selected
		elif not graph.eventGraph.testConstraint(graph.segmentThreads[other.segment], graph.segmentThreads[side.segment]):
			selected.remove(other)
	return selected

def oldest_sides(graph):
	return reduce(lambda X, Y: select_oldest_sides(X, Y, graph), graph.sides(), [])

###############################################
## Case 1
## Adding G-reducible necessary bridge labels
###############################################

def listCase1(graph):
	pass

def applyCase1(args):
	pass

###############################################
## Case 2
## Adding a G-reducible oldest junction label whose most ancestral 
## labelled descendants are junctions or G-irreducible
###############################################

def listCase2(graph):
	return Extension(data, applyCase2)

def applyCase2(args):
	pass

###############################################
## Case 3
## Adding G-reducible necessary bridge bonds
###############################################

def listCase3(graph):
	pass

def applyCase3(args):
	pass

###############################################
## Case 4
## Adding a G-reducible oldest bond whose most ancestral 
## labelled descendants are junctions or G-irreducible
###############################################

def listCase4(graph):
	pass

def applyCase4(args):
	pass

###############################################
## Case 5
## The addition of a G-reducible free tailed (non-minimal) branch
###############################################

def listCase5(graph):
	pass

def applyCase5(args):
	pass

detectors = [listCase1, listCase2, listCase3, listCase4, listCase5]
