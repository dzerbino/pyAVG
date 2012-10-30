#!/usr/bin/env python

import segment
import copy

class Traversal(object):
	def __init__(self, seg, orientation):
		assert isinstance(seg, segment.Segment)
		self.segment = seg
		self.orientation = bool(orientation)

	def __copy__(self):
		return Traversal(copy.copy(self.segment), self.orientation)

	def dot(self):
		# Could be nice to change display w/ orientation
		return self.segment.dot()

	def opposite(self):
		return Traversal(self.segment, not self.orientation)

	def connect(self, other):
		if self.orientation and other.orientation:
			self.segment.right.createBond(other.segment.left)
		elif self.orientation and not other.orientation:
			self.segment.right.createBond(other.segment.right)
		elif not self.orientation and other.orientation:
			self.segment.left.createBond(other.segment.left)
		else:
			self.segment.left.createBond(other.segment.right)

	def isConnected(self, other):
		if self.orientation and other.orientation:
			return self.segment.right is other.segment.left.bond
		elif self.orientation and not other.orientation:
			return self.segment.right is other.segment.right.bond
		elif not self.orientation and other.orientation:
			return self.segment.left is other.segment.left.bond
		else:
			return self.segment.left is other.segment.right.bond

