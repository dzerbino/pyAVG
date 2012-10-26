#!/usr/bin/env python

class DNAHistoryGraph(object):
	def __init__(self, segments):
		self.segments = list(segments)
	
	def isAVG(self):
		return not any(X.isAmbiguous() for X in self.segments)

	def sequenceComplexity(self, lowerBound=True):
		return sum(X.sequenceComplexity(lowerBound) for X in self.segments)

	def rearrangementComplexity(self, lowerBound=True):
		return sum(X.rearrangementComplexity(lowerBound) for X in self.modules())
