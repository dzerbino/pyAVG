#!/usr/bin/env python

class Label(object):
	def __init__(self, segment, sequence=None):
		self.segment = segment
		self.sequence = str(sequence)[0]

	def __str__(self):
		return self.sequence

	def complement(self):
		return { "A":"T", "T":'A', "C":'G', "G":"C"}[self.sequence]

	def __cmp__(self, other):
		if other is None:
			return 1
		else:
			return cmp(self.sequence, other.sequence)

	def __hash__(self):
		return self.sequence.__hash__()
