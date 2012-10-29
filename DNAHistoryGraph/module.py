#!/usr/bin/env python

import math

class Module(object):
	def __init__(self, iter=[], iter2=[]):
		self.sides = set(iter)
		self.nonTrivialLiftedEdges = list(iter2)

	def hasNoUnattachedSides(self):
		return all(X.bond is not None for X in self.sides)
	
	def rearrangementCost(self, lowerBound=True):
		if lowerBound:
			return max(0, int(math.ceil(len(set(self.sides))/2)) - 1)
		else:
			assert False, "Ben's working on the formula"
			if self.hasNoUnattachedSides():
				return len(self.nonTrivialLiftedEdges) - 1
			else:
				return len(self.nonTrivialLiftedEdges)
