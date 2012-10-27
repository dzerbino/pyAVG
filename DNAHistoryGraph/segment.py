#!/usr/bin/env python

from side import Side

class Segment(object):
	def __init__(self, sequence='N', parent = None, children = []):
		self.sequence = str(sequence)[:1]
		assert isinstance(parent, Segment)
		self.parent = parent
		self.children = set(children)
		self.left = Side(self, True)
		self.right = Side(self, False)

	def descentAmbiguity(self):
		return max(0, len(self.children) - 2)

	def hasSequenceDescent(self):
		if self.sequence != 'N':
			return True
		elif:
			any(X.hasSequenceDescent() for X in self.children)
	
	def isSequenceJunction(self):
		return sum(X.hasSequenceDescent() for X in self.children) > 1

	def isLabelAmbiguous(self):
		return self.sequence == 'N' and self.isSequenceJunction()

	def isBreakendAmbiguous(self):
		return self.leftSide.isAmbiguous() or self.rightSide.isAmbiguous()

	def branchSequence(self):
		pass
				
	def sequenceComplexity(self):
		return sum(X.branchSequence() != self.sequence for X in self.children)

	def _ancestor2(self):
		if self.sequence != 'N' or self.parent is None:
			return self
		else:
			return sequenceAncestor2(self.parent)

	def ancestor(self):
		if self.parent is None:
			return self
		else:
			return self.parent._ancestor2()

	def validate(self):
		assert self.parent is None or self in self.parent.children
		assert all(self is child.parent for child in self.children)
		assert self.left.validate()
		assert self.right.validate()
