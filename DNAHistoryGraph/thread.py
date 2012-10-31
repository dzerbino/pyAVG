#!/usr/bin/env python

import segment
import traversal
import copy

class Thread(object):
	""" Walk through DNA history graph, i.e. sequence of oriented segments """
	################################
	## Basics
	################################
	def __init__(self, iter=[]):
		""" Creates a thread from a sequence of traversals, connects them up as described, then expands if possible """
		self.traversals = list(iter)
		self._connect()
		self._expand()

	def __getitem__(self, key):
		return self.traversals[key]

	def __len__(self):
		return len(self.traversals)

	def __copy__(self):
		thread = Thread(copy.copy(X) for X in self)
		thread._connect()
		if self[-1].isConnected(self[0]):
			thread[-1].connect(thread[0])
		return thread

	def append(self, element):
		self.traversals.append(element)

	def segments(self):
		return [X.segment for X in self.traversals]

	def _connect(self):
		for traversalA, traversalB in zip(self[:-1], self[1:]):
			traversalA.connect(traversalB)

	################################
	## Output
	################################

	def dot(self):
		ranks = " ".join(["{rank = same;"] + [str(id(X.segment)) for X in self.traversals] + ["}"])
		data = "\n".join(X.dot() for X in self.traversals)
		return "\n".join([ranks, data])

	################################
	## Operations
	################################

	def isCycle(self):
		return self[-1].isConnected(self[0])

	def childThread(self):
		""" Produces thread of new children nodes of the thread nodes """
		T = copy.copy(self)
		for traversalA, traversalB in zip(self, T):
			traversalA.createBranch(traversalB)
		return T

	def _expandRight(self):
		""" Expand a thread to its right, following as far as possible, stopping when a cycle is created """
		if not self.isCycle():
			next = self[-1].next()
			if next is not None:
				self.traversals.append(next)
				self._expandRight()
			
	def _expandLeft(self):
		""" Expand a thread to its left, following bonds as far as possible, stopping when a cycle is created """
		if not self.isCycle():
			previous = self[0].previous()
			if previous is not None:
				self.traversals = [previous] + self.traversals
				self._expandLeft()

	def _expand(self):
		self._expandRight()
		self._expandLeft()

	def sequence(self):
		""" Returns sequence assigned to thread """
		return "".join(lambda X: X.sequence() for X in self)

	################################
	## Validation
	################################
	def validate(self):
		if len(self) > 0:
			copy = self[0].segment.thread()
			assert len(self) == len(copy), "%s\t%s" % (len(self), len(copy))
			selfSegments = self.segments()
			copySegments = copy.segments()
			assert all(X[0] == X[1] for X in zip(sorted(selfSegments), sorted(copySegments)))

class CircularThread(Thread):
	""" Circular thread """
	def __init__(self, iter):
		elems = list(iter)
		if len(elems) > 0:
			elems[-1].connect(elems[0])
		super(CircularThread, self).__init__(elems)

class SequenceThread(Thread):
	""" Thread created from a string sequence """
	def __init__(self, sequence):
		super(SequenceThread, self).__init__(traversal.Traversal(segment.Segment(sequence = X), True) for X in str(sequence))

class CircularSequenceThread(CircularThread):
	""" Circular thread created from a string sequence """
	def __init__(self, sequence):
		super(CircularSequenceThread, self).__init__(traversal.Traversal(segment.Segment(sequence = X), True) for X in str(sequence))
	
