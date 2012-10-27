#!/usr/bin/env python

from segment import Segment

class Side(object):
	def __init__(self, segment, left):
		assert isinstance(segment, Segment)
		self.segment = segment
		self.left = bool(left)

	def opposite(self):
		if self.left:
			return self.segment.right
		else:
			return self.segment.left

	def parent(self):
		if self.segment.parent is not None:
			if self.left:
				return self.segment.parent.left
			else:
				return self.segment.parent.right
		else:
			return None

	def children(self):
		if self.left:
			return X.left for X in self.segment.children
		else:
			return X.right for X in self.segment.children

	def _isAttached(self):
		if self.left:
			return self.segment.leftBond is not None
		else:
			return self.segment.rightBond is not None

	def _ancestor2(self):
		if self._isAttached() or self.parent() is None:
			return self
		else:
			return self.parent._ancestor2()

	def ancestor(self):
		if self.parent() is None:
			return self.parent()._ancestor2()
		else:
			return self

	def hasSideBreakendDescent(self, left):
		if left and self.leftBond is not None:
			return True
		elif not left and self.rightBond is not None:
			return False
		else:
			return any(X.hasSideBreakendDescent(left) for X in self.children)

	def isAmbiguous(self):
		if self._isAttached:
			return False
		else:
			return sum(X.hasSideBreakendDescent() for X in self.children()) > 1

	def createBond(self, other):
		self.bond = other
		other.bond = self

	def deleteBond(self):
		if self.bond is not None:
			self.bond.bond = None
			self.bond = None

	def validate(self):
		assert self.bond is None or self.bond.bond is self
