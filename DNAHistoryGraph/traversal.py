#!/usr/bin/env python

import segment
import copy

class Traversal(object):
	""" 
	A segment as seen from a given traversal orientation. 

	If the orientation is direct (i.e. True), the traversal starts on the left side of the segment and ends on the right.

	If the orientation is opposite (i.e. False), the traversal starts on the right side of the segment and ends on the left.
	"""
	######################################
	## Basics
	######################################
	def __init__(self, seg, orientation):
		assert isinstance(seg, segment.Segment)
		self.segment = seg
		self.orientation = bool(orientation)

	def __copy__(self):
		return Traversal(copy.copy(self.segment), self.orientation)

	def __eq__(self, other):
		return self.segment == other.segment and self.orientation == other.orientation

	######################################
	## Output
	######################################
	def dot(self):
		# Could be nice to change display w/ orientation
		return self.segment.dot()

	######################################
	## Convenience functions
	######################################
	def opposite(self):
		""" The traveral of the same segment but in the opposite direction """
		return Traversal(self.segment, not self.orientation)

	def start(self):
		""" Start side of the traversed segment """
		if self.orientation:
			return self.segment.left
		else:
			return self.segment.right

	def previous(self):
		""" Returns traversal of the segment prior to the current segment traversal """
		bond = self.start().bond
		if bond is None:
			return None
		else:
			return Traversal(bond.segment, not bond.left)

	def end(self):
		""" End side of the traversed segment """
		return self.start().opposite

	def next(self):
		""" Returns traversal of the segment following the current segment traversal """
		bond = self.end().bond
		if bond is None:
			return None
		else:
			return Traversal(bond.segment, bond.left)

	def connect(self, other):
		""" Connect the end of traversal to the beginning of another traversal. """
		self.end().createBond(other.start())

	def isConnected(self, other):
		""" Check whether the end of a traversal is connected to the beginning of another traversal. """
		return self.end().bond is other.start()

	def createBranch(self, other):
		""" Create branch from segment to segment of other traversal """
		self.segment.children.add(other.segment)
		other.segment.parent = self.segment

	def sequence(self):
		""" Returns sequence of segment, reverse complemented if the orientation is negative """
		if self.segment.label is None:
			return 'N'
		elif self.orientation:
			return str(self.segment.label)
		elif self.segment.label is None:
			return self.segment.label.complement()
