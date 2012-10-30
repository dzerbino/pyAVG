#!/usr/bin/env python

import segment
import copy

def _invert(segment):
	return (segment[0], not segment[1])

def _connectPair(segmentA, segmentB):
	if segmentA[1] and segmentB[1]:
		segmentA[0].right.createBond(segmentB[0].left)
	elif segmentA[1] and not segmentB[1]:
		segmentA[0].right.createBond(segmentB[0].right)
	elif not segmentA[1] and segmentB[1]:
		segmentA[0].left.createBond(segmentB[0].left)
	else:
		segmentA[0].left.createBond(segmentB[0].right)

def _pairIsConnected(segmentA, segmentB):
	if segmentA[1] and segmentB[1]:
		return segmentA[0].right is segmentB[0].left.bond
	elif segmentA[1] and not segmentB[1]:
		return segmentA[0].right is segmentB[0].right.bond
	elif not segmentA[1] and segmentB[1]:
		return segmentA[0].left is segmentB[0].left.bond
	else:
		return segmentA[0].left is segmentB[0].right.bond

def redirect(segmentA, invert, segmentB, invert2):
	if invert:
		segmentA = _invert(segmentA)
	if invert2:
		segmentB = _invert(segmentB)
	_connectPair(segmentA, segmentB)

class Thread(object):
	""" Walk through DNA history graph, i.e. sequence of oriented segments """
	def __init__(self, iter=[]):
		self.elements = list(iter)
		assert all(isinstance(X[0], segment.Segment) for X in self)

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

	def segments(self):
		return [X[0] for X in self.elements]

	def dot(self):
		ranks = " ".join(["{rank = same;"] + [str(id(X[0])) for X in self.elements] + ["}"])
		data = "\n".join(X[0].dot() for X in self.elements)
		return "\n".join([ranks, data])

	def childThread(self):
		T = copy.copy(self)
		for segmentA, segmentB in zip(self, T):
			segmentA[0].children.add(segmentB[0])
			segmentB[0].parent = segmentA[0]
		return T

	def expandRight(self):
		tail, orientation = self[-1]
		if orientation:
			bond = tail.right.bond
			if bond is None or bond.segment in self.segments():
				return
			else:
				self.elements.append((bond.segment, bond.left))
				return self.expandRight()
		else:
			bond = tail.left.bond
			if bond is None or bond.segment in self.segments():
				return
			else:
				self.elements.append((bond.segment, bond.left))
				return self.expandRight()
			
	def expandLeft(self):
		head, orientation = self[0]
		if orientation:
			bond = head.left.bond
			if bond is None or bond.segment in self.segments():
				return
			else:
				self.elements = [(bond.segment, not bond.left)] + self.elements
				return self.expandLeft()
		else:
			bond = head.right.bond
			if bond is None or bond.segment in self.segments():
				return
			else:
				self.elements = [(bond.segment, not bond.left)] + self.elements
				return self.expandLeft()

	def redirect(self, pos1, invert, pos2, invert2):
		redirect(self[pos1 % len(self)], invert, self[pos2 % len(self)], invert2)

	def validate(self):
		if len(self) > 0:
			copy = self[0][0].thread()
			assert all(X in self.segments() for X in copy.segments())
			assert all(X in copy.segments() for X in self.segments())

class SequenceThread(Thread):
	def __init__(self, sequence):
		super(SequenceThread, self).__init__((segment.Segment(sequence = X), True) for X in str(sequence))
		self._connect()

class CircularSequenceThread(SequenceThread):
	def __init__(self, sequence):
		super(CircularSequenceThread, self).__init__(sequence)
		self[0][0].left.createBond(self[-1][0].right)
