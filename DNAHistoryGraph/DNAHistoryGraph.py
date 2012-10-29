#!/usr/bin/env python

from segment import Segment
from partialOrderSet import PartialOrderSet

class DNAHistoryGraph(object):
	""" DNA History graph """

	##################################
	## Basics
	##################################
	def __init__(self, segments=[]):
		self.segments = list(segments)
		assert all(isinstance(X, Segment) for X in self.segments)
		print len(self.segments), 'INIT'
		self.eventGraph, self.segmentThreads = self.threads()
		self.timeEventGraph()

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
			thread = sideA.segment.expandThread(Thread())

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
			thread = sideA.segment.expandThread(Thread())
			if sideB not in thread:
				oldThread = self.segmentThreads[sideA]
				thread2 = sideB.segment.expandThread(Thread())

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
		print 'HIST'
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
		assert self.eventGraph.validate()
		return True
