#!/usr/bin/env python

class Label(object):
	def __init__(self, segment, sequence=None):
		self.segment = segment
		self.sequence = str(sequence)[0]

	def __str__(self):
		return self.sequence

	def complement(self):
		if self.sequence is 'A':
			return 'T'
		else:
			return 'A'

	def __cmp__(self, other):
		if other is None:
			return 1
		else:
			return cmp(self.sequence, other.sequence)

	def __hash__(self):
		return id(self)
