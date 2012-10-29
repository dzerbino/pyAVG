#!/usr/bin/env python

from segment import Segment

class Thread(object):
	def __init__(self, iter = None):
		self.segments = list(iter)
		self._connectNodes()

	def _connectNodes(self):
		for segmentA, segmentB in zip(self.segments[:-1], self.segments[1:]):
			segmentA.rightBond = segmentB.leftSide()
			segmentB.leftBond = segmentA.rightSide()

class SequenceThread(Thread):
	def __init__(self, sequence):
		super(SequenceThread, self).__init__(Segment(sequence = X) for X in str(sequence))

class CircularThread(Thread):
	def __init__(self, iter = None):
		super(CircularThread, self).__init__(iter)
		self.segments[0].leftBond = self.segments[-1].rightSide()
		self.segments[-1].rightBond = self.segments[0].leftSide()

class CircularSequenceThread(CircularThread):
	def __init__(self, sequence):
		super(CircularSequenceThread, self).__init__(Segment(sequence = X) for X in str(sequence))
