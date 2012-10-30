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
		self.traversals = list(iter)
		self._connect()
		assert all(isinstance(X, traversal.Traversal) for X in self)

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

	def childThread(self):
		T = copy.copy(self)
		for traversalA, traversalB in zip(self, T):
			traversalA.segment.children.add(traversalB.segment)
			traversalB.segment.parent = traversalA.segment
		return T

	def expandRight(self):
		tail = self[-1]
		if tail.orientation:
			bond = tail.segment.right.bond
			if bond is None or bond.segment in self.segments():
				return
			else:
				self.traversals.append(traversal.Traversal(bond.segment, bond.left))
				return self.expandRight()
		else:
			bond = tail.segment.left.bond
			if bond is None or bond.segment in self.segments():
				return
			else:
				self.traversals.append(traversal.Traversal(bond.segment, bond.left))
				return self.expandRight()
			
	def expandLeft(self):
		head = self[0]
		if head.orientation:
			bond = head.segment.left.bond
			if bond is None or bond.segment in self.segments():
				return
			else:
				self.traversals = [traversal.Traversal(bond.segment, not bond.left)] + self.traversals
				return self.expandLeft()
		else:
			bond = head.segment.right.bond
			if bond is None or bond.segment in self.segments():
				return
			else:
				self.traversals = [traversal.Traversal(bond.segment, not bond.left)] + self.traversals
				return self.expandLeft()

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
		super(CircularThread, self).__init__(iter)
		if len(self) > 0:
			self[-1].connect(self[0])

class SequenceThread(Thread):
	""" Thread created from a string sequence """
	def __init__(self, sequence):
		super(SequenceThread, self).__init__(traversal.Traversal(segment.Segment(sequence = X), True) for X in str(sequence))

class CircularSequenceThread(CircularThread):
	""" Circular thread created from a string sequence """
	def __init__(self, sequence):
		super(CircularSequenceThread, self).__init__(traversal.Traversal(segment.Segment(sequence = X), True) for X in str(sequence))
	
