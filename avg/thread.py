#!/usr/bin/env python

class Thread(object):
	def __init__(self, iter = None):
		self.segments = list(iter)
		self._connectNodes()

	def _connectNodes(self):
		for segmentA, segmentB in zip(self.segments[:-1], self.segments[1:]):
			segmentA.rightBond = Side(segmentB, True)
			segmentB.leftBond = Side(segmentA, False)


class SequenceThread(Thread):
	def __init__(self, sequence):
		super(SequenceThread, self).__init__(Segment(sequence = X) for X in str(sequence))
