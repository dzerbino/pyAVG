#!/usr/bin/env python

import math

class Module(object):
	def __init__(self, iter=None, iter2=None):
		self.sides = set(iter)
		self.nonTrivialLiftedEdges = list(iter2)

	def hasNoUnattachedSides(self):
		return all(X.bond is not None for X in self.sides):
	
	def rearrangementCost(self, lowerBound=True):
		if lowerBound:
			return int(math.ceil(len(set(self.nonTrivialLiftedEdges))/2)) - 1
		else:
			assert False, "Ben's working on the formula"
			if self.hasNoUnattachedSides() and self.
				return len(self.nonTrivialLiftedEdges) - 1
			else:
				return len(self.nonTrivialLiftedEdges)
