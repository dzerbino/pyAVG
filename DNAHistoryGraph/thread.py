#!/usr/bin/env python

import segment
import copy

def _connectPair(segmentA, segmentB):
	if segmentA[1] and segmentB[1]:
		segmentA[0].right.createBond(segmentB[0].left)
	elif segmentA[1] and not segmentB[1]:
		segmentA[0].right.createBond(segmentB[0].right)
	elif not segmentA[1] and segmentB[1]:
		segmentA[0].left.createBond(segmentB[0].right)
	else:
		segmentA[0].left.createBond(segmentB[0].left)

def _pairIsConnected(segmentA, segmentB):
	if segmentA[1] and segmentB[1]:
		return segmentA[0].right is segmentB[0].left
	elif segmentA[1] and not segmentB[1]:
		return segmentA[0].right is segmentB[0].right
	elif not segmentA[1] and segmentB[1]:
		return segmentA[0].left is segmentB[0].right
	else:
		return segmentA[0].left is segmentB[0].left

class Thread(object):
	""" Walk through DNA history graph, i.e. sequence of oriented segments """
	def __init__(self, iter=[]):
		self.elements = list(iter)

	def __getitem__(self, key):
		return self.elements[key]

	def __len__(self):
		return len(self.elements)

	def _connect(self):
		for segmentA, segmentB in zip(self[:-1], self[1:]):
			_connectPair(segmentA, segmentB)

	def __copy__(self):
		thread = Thread((copy.copy(X[0]), X[1]) for X in self)
		thread._connect()
		if _pairIsConnected(self[-1], self[0]):
			_connectPair(thread[-1], thread[0])
		return thread

	def append(self, value, orientation):
		self.elements.append(value, orientation)

	def dot(self):
		ranks = " ".join(["{rank = same;"] + [str(id(X)) for X in self.elements] + ["}"])
		data = "\n".join(X.dot() for X in self.elements)
		return "\n".join([ranks, data])

	def expandRight(self):
		assert len(self) > 0
		tail, orientation = self[-1]
		if orientation:
			bond = tail.right.bond
			if bond is None or bond.segment is self.elements[0][0]:
				return
			else:
				self.elements.append((bond.segment, bond.orientation))
				return self.expandRight()
		else:
			bond = tail.left.bond
			if bond is None or bond.segment is self.elements[0][0]:
				return
			else:
				self.elements.append((bond.segment, bond.orientation))
				return self.expandRight()
			
	def expandLeft(self):
		assert len(self) > 0
		tail, orientation = self[-1]
		if orientation:
			bond = tail.left.bond
			if bond is None:
				return
			else:
				self.elements.append((bond.segment, not bond.orientation))
				return self.expandLeft()
		else:
			bond = tail.right.bond
			if bond is None:
				return
			else:
				self.elements.append((bond.segment, not bond.orientation))
				return self.expandLeft()

class SequenceThread(Thread):
	def __init__(self, sequence):
		super(SequenceThread, self).__init__((segment.Segment(sequence = X), True) for X in str(sequence))
		self._connect()

class CircularSequenceThread(SequenceThread):
	def __init__(self, sequence):
		super(CircularSequenceThread, self).__init__(sequence)
		self[0][0].left.createBond(self[-1][0].right)
