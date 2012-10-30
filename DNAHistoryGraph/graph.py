#!/usr/bin/env python

import copy

from segment import Segment
from pyAVG.utils.partialOrderSet import PartialOrderSet

class DNAHistoryGraph(object):
	""" DNA History graph """

	##################################
	## Basics
	##################################
	def __init__(self, segments=[]):
		self.segments = list(segments)
		self.eventGraph, self.segmentThreads = self.threads()
		self.timeEventGraph()

	def __copy__(self):
		duplicates = dict(map(lambda X: (X, copy.copy(X)), self.segments))
		for segment in self.segments:
			if segment.left.bond is not None:
				if segment.left.bond.left:
					duplicates[segment].left.createBond(duplicates[segment.left.bond.segment].left)
				else:
					duplicates[segment].left.createBond(duplicates[segment.left.bond.segment].right)
			if segment.right.bond is not None:
				if segment.right.bond.left:
					duplicates[segment].right.createBond(duplicates[segment.right.bond.segment].left)
				else:
					duplicates[segment].right.createBond(duplicates[segment.right.bond.segment].right)

			if segment.parent is not None:
				duplicates[segment].parent = duplicates[segment.parent]
			duplicates[segment].children = set(duplicates[X] for X in segment.children)
		return DNAHistoryGraph(duplicates.values())

	##################################
	## Online acyclicity verification
	##################################
	def threads(self):
		""" Computes tuples (PartialOrderSet X, dict Y) which contains X) graph threads (no ordering) and Y) segment to thread mapping """
		return reduce(lambda X, Y: Y.threads(X), self.segments, (PartialOrderSet(), dict()))

	def timeEventGraph(self):
		""" Adds timing constraints to unordered set of threads """	
		assert all(X.parent in self.segments for X in self.segments if X.parent is not None)
		assert all(X.parent in self.segmentThreads for X in self.segments if X.parent is not None)
		assert all(X in self.segmentThreads for X in self.segments)
		for segment in self.segments:
			if segment.parent is not None:
				self.eventGraph.addConstraint(self.segmentThreads[segment.parent], self.segmentThreads[segment])

	def createBond(self, sideA, sideB):
		""" Creates bond between two sides and throws RuntimeError if cycle created """
		if sideA.bond is not None and sideA.bond is not sideB:
			self.deleteBond(sideA)
		sideA.createBond(sideB)
		if sideB.segmentThread[sideB] is not self.segmentThreads[sideA]:
			oldThread = self.segmentThreads[sideA]
			oldThread2 = self.segmentThreads[sideB]
			thread = sideA.segment.thread()

			# Updating self.segmentThreads
			for segment in thread:
				self.segmentThreads[segment] = thread

			# Updating self.eventGraph
			self.eventGraph.pop(oldThread)
			self.eventGraph.pop(oldThread2)
			self.eventGraph.addElement(thread)
			for segment in thread:
				self.eventGraph.addConstraint(self.segmentThreads[segment.parent], thread)
				for child in segment.children:
					self.eventGraph.addConstraint(thread, self.segmentThreads[child])
		
	def deleteBond(self, sideA):
		""" Deletes bond between two sides and updates event graph """
		sideB = sideA.bond
		sideA.deleteBond()
		if sideB is not None:
			thread = sideA.segment.thread()
			if sideB not in thread:
				oldThread = self.segmentThreads[sideA]
				thread2 = sideB.segment.thread()

				# Updating self.segmentThreads
				for segment in thread:
					self.segmentThreads[segment] = thread
				for segment in thread2:
					self.segmentThreads[segment] = thread2

				# Updating self.eventGraph
				self.eventGraph.pop(oldThread)
				self.eventGraph.addElement(thread)
				self.eventGraph.addElement(thread2)
				for segment in thread:
					self.eventGraph.addConstraint(self.segmentThreads[segment.parent], thread)
					for child in segment.children:
						self.eventGraph.addConstraint(thread, self.segmentThreads[child])
				for segment in thread2:
					self.eventGraph.addConstraint(self.segmentThreads[segment.parent], thread2)
					for child in segment.children:
						self.eventGraph.addConstraint(thread2, self.segmentThreads[child])


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
		return self.coalescenceAmbiguity() + self.substitutionAmbiguity() + self.rearrangementAmbiguity()

	def isAVG(self):
		return self.ambiguity() == 0

	##################################
	## Cost
	##################################
	def substitutionCost(self, lowerBound=True):
		return sum(X.substitutionCost(lowerBound) for X in self.segments)

	def modules(self):
		return reduce(lambda X, Y: Y.modules(X), self.segments, (list(), set()))[0]

	def rearrangementCost(self, lowerBound=True):
		return sum(X.rearrangementCost(lowerBound) for X in self.modules())

	##################################
	## Output
	##################################
	def dot(self):
		return "\n".join(["digraph G {"] + [X.dot() for X in self.eventGraph] + ["}"])

	##################################
	## Validation
	##################################
	def validate(self):
		assert all(X.validate() for X in self.segments)
		assert all(X.parent in self.segments for X in self.segments if X.parent is not None)
		assert all(Y in self.segments for X in self.segments for Y in X.children)
		assert all(X.left.bond.segment in self.segments for X in self.segments if X.left.bond is not None)
		assert all(X.right.bond.segment in self.segments for X in self.segments if X.right.bond is not None)
		assert self.eventGraph.validate()
		return True
