#!/usr/bin/env python

import math

class Module(object):
	def __init__(self, iter=[], iter2=[]):
		self.sides = set(iter)
		self.nonTrivialLiftedEdges = list(iter2)

	def _cycleDiscount(self):
		if True:
			return 1
		else:
			return 0
	
	def rearrangementCost(self, lowerBound=True):
		""" Returns lower bound (default) or upper bound to the rearrangement cost of a module """
		if lowerBound:
			return max(0, int(math.ceil(len(set(self.sides))* 0.5)) - 1)
		else:
			return len(self.nonTrivialLiftedEdges) - self._cycleDiscount()
