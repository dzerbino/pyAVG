#!/usr/bin/env python

from segment import Segment
from partialOrderSet import PartialOrderSet

def createEventGraph(segments):

class DNAHistoryGraph(object):
	""" DNA History graph """

	##################################
	## Basics
	##################################
	def __init__(self, segments):
		assert all(isinstance(X, Segment) for X in segments)
		self.segments = list(segments)
		self.threads, self.segmentThreads = self.threads()
		self.eventGraph = self.createEventGraph()

	##################################
	## Online acyclicity verification
	##################################
	def threads(self):
		return reduce(lambda X, Y: Y.threads(X), self.segments, (list(), set()))[0]

	def createEventGraph(self):
		eventGraph = PartialOrderSet(self.threads)
		for segment in self.segments:
			assert self.addConstraint(self.segmentThreads[segment.parent], self.segmentThreads[segment])
		return eventGraph

	##################################
	## Ambiguity
	##################################
	def substitutionAmbiguity(self):
		return sum(segment.substitutionAmbiguity() for segment in self.segments)

	def coalescenceAmbiguity(self):
		return sum(segment.coalescenceAmbiguity() for segment in self.segments)

	def rearrangementAmbiguity(self):
		return sum(segment.rearrangementAmbiguity() for segment in self.segments)
	
	def ambiguity(self):
		return self.descentAmbiguity() + self.labelAmbiguity() + self.rearrangementAmbiguity()

	def isAVG(self):
		return self.ambiguity() == 0

	##################################
	## Cost
	##################################
	def substitutionCost(self, lowerBound=True):
		return sum(X.substitutionCost(lowerBound) for X in self.segments)

	def modules():
		return reduce(lambda X, Y: Y.modules(X), self.segments, (list(), set()))[0]

	def rearrangementCost(self, lowerBound=True):
		return sum(X.rearrangementCost(lowerBound) for X in self.modules())

	##################################
	## Validation
	##################################
	def validate(self):
		assert all(X.validate() for X in self.segments)
		assert self.eventGraph.validate()
