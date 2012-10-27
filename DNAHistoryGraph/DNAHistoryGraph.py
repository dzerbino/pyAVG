#!/usr/bin/env python

from segment import Segment

def createEventGraph(segments):

class DNAHistoryGraph(object):
	def __init__(self, segments):
		assert all(isinstance(X, Segment) for X in segments)
		self.segments = list(segments)
		self.eventGraph = createEventGraph(segments)

	def labelAmbiguity(self):
		return sum(segment.isLabelAmbiguous() for segment in self.segments)

	def descentAmbiguity(self):
		return 2 * sum(segment.descentAmbiguity() for segment in self.segments)
	
	def ambiguity(self):
		return self.descentAmbiguity() + self.labelAmbiguity()

	def isAVG(self):
		return self.ambiguity() == 0

	def sequenceComplexity(self, lowerBound=True):
		return sum(X.sequenceComplexity(lowerBound) for X in self.segments)

	def rearrangementComplexity(self, lowerBound=True):
		return sum(X.rearrangementComplexity(lowerBound) for X in self.modules())

	def validate(self):
		assert all(X.validate() for X in self.segments)
		assert self.eventGraph.validate()
