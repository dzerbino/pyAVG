#!/usr/bin/env python

import segment

class Thread(object):
	def __init__(self, iter=[]):
		self.elements = list(iter)

	def __getitem__(self, key):
		return self.elements[key]

	def __len__(self):
		return len(self.elements)

	def append(self, value):
		self.elements.append(value)

	def dot(self):
		ranks = " ".join(["{rank = same;"] + [str(id(X)) for X in self.elements] + ["}"])
		data = "\n".join(X.dot() for X in self.elements)
		return "\n".join([ranks, data])

class SequenceThread(Thread):
	def __init__(self, sequence):
		super(SequenceThread, self).__init__(segment.Segment(sequence = X) for X in str(sequence))
		self._connectNodes()

	def _connectNodes(self):
		for segmentA, segmentB in zip(self[:-1], self[1:]):
			segmentA.right.createBond(segmentB.left)

class CircularThread(SequenceThread):
	def __init__(self, iter = None):
		super(CircularThread, self).__init__(iter)
		self[0].left.createBond(self[-1].right)

class CircularSequenceThread(CircularThread):
	def __init__(self, sequence):
		super(CircularSequenceThread, self).__init__(sequence)
