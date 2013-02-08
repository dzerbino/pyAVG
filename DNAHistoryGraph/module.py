#!/usr/bin/env python

import math
import operator

class Module(object):
	def __init__(self, side):
		self.sides = set()
		self.freeRootNumber = 0
		def expandModule(side):
			if side != None and side not in self.sides:
				assert side.isModuleMaterial()
				self.sides.add(side)
				expandModule(side.bond)
				for descendant in side.nonTrivialLiftedBonds():
					if descendant.bond != descendant.bond.ancestor():
						expandModule(descendant.bond.ancestor())
					else:
						self.freeRootNumber += 1
		expandModule(side)
	
	def lowerBoundRearrangementCost(self):
		#print "calc", self.freeRootNumber, len(self.sides), math.ceil((self.freeRootNumber + len(self.sides))/2.0) - 1
		return math.ceil((self.freeRootNumber + len(self.sides))/2.0) - 1
	
	def upperBoundRearrangementCost(self):
		return (sum([ len(x.nonTrivialLiftedBonds()) for x in self.sides ]) + self.freeRootNumber)/2 - reduce(operator.mul, [ len(x.nonTrivialLiftedBonds()) == 1 for x in self.sides ], 1)

	def _liftedEdgeDot(self, liftedEdge):
		sides = list(liftedEdge)
		return "%s -> %s [style=dashed]" % (id(sides[0].segment), id(sides[1].segment))

	def _liftedEdgesDot(self):
		return map(self._liftedEdgeDot, self.nonTrivialLiftedEdges)

	def dot(self):
		return "\n".join(["node [color=blue]"] + map(lambda X: str(id(X.segment)), self.sides) + self._liftedEdgesDot() + ["node [color=black]"])

	def isSimple(self):
		return reduce(operator.mul, [ len(x.nonTrivialLiftedBonds()) <= 1 for x in self.sides ], 1) == 1

	def validate(self, graph):
		assert all(X.bond is not None or X.parent() is None for X in self.sides)
		assert self.lowerBoundRearrangementCost() <= self.upperBoundRearrangementCost(), "\n".join([graph.dot(), self.dot(), str(self.lowerBoundRearrangementCost()), str(self.upperBoundRearrangementCost())])
		return True
