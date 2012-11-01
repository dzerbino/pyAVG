#!/usr/bin/env python

import math

class Module(object):
	def __init__(self, iter=[], iter2=[]):
		self.sides = set(iter)
		self.nonTrivialLiftedEdges = list(iter2)

	def _hasNoUnattachedSides(self):
		return not any(X.bond is None for X in self.sides)

	def _hasALiftedEdgeForEveryTwoSides(self):
		return len(self.sides) == 2 * len(self.nonTrivialLiftedEdges)

	def _cycleDiscount(self):
		if self._hasNoUnattachedSides() and self._hasALiftedEdgeForEveryTwoSides():
			return 1
		else:
			return 0
	
	def rearrangementCost(self, lowerBound=True):
		""" Returns lower bound (default) or upper bound to the rearrangement cost of a module """
		if lowerBound:
			return max(0, int(math.ceil(len(set(self.sides))* 0.5)) - 1)
		else:
			return len(self.nonTrivialLiftedEdges) - self._cycleDiscount()

	def _liftedEdgeDot(self, liftedEdge):
		sides = list(liftedEdge)
		return "%s -> %s [style=dashed]" % (id(sides[0].segment), id(sides[1].segment))

	def _liftedEdgesDot(self):
		return map(self._liftedEdgeDot, self.nonTrivialLiftedEdges)

	def dot(self):
		return "\n".join(["node [color=blue]"] + map(lambda X: str(id(X.segment)), self.sides) + self._liftedEdgesDot() + ["node [color=black]"])

	def isSimple(self):
		# Test whether there are no duplicates in list of lifted edge nodes
		sides = sum([list(X) for X in self.nonTrivialLiftedEdges], [])
		return len(sides) == len(set(sides))

	def validate(self, graph):
		assert all(X.bond is not None or X.parent() is None for X in self.sides)
		assert self.rearrangementCost(lowerBound=True) <= self.rearrangementCost(lowerBound=False), "\n".join([graph.dot(), self.dot(), str(self.rearrangementCost(lowerBound=True)), str(self.rearrangementCost(lowerBound=False))])
		return True
