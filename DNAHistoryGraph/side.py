#!/usr/bin/env python

class Side(object):
	def __init__(self, segment, left):
		self.segment = segment
		self.left = bool(left)

	def parent(self):
		if self.segment.parent is not None:
			return Side(self.segment.parent, self.left)
		else:
			return None

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
		return self.parent()._ancestor2()

	def hasSideBreakendDescent(self, left):
		if left and self.leftBond is not None:
			return True
		elif not left and self.rightBond is not None:
			return False
		else:
			return any(X.hasSideBreakendDescent(left) for X in self.children)

	def isSideBreakendAmbiguous(self, left):
		if left and self.leftBond is not None:
			return False
		elif not left and self.rightBond is not None:
			return False
		else:
			return sum(X.hasSideBreakendDescent(left) for X in self.children) > 1
