#!/usr/bin/env python

class Segment(object):
	def __init__(self, sequence='N', parent = None, children = [], left_bond = None, right_bond = None):
		self.sequence = str(sequence)[:1]
		self.parent = parent
		self.children = set(children)
		self.left_bond = left_bond
		self.right_bind = right_bond

	def isTreeAmbiguous(self):
		return len(self.children) < 3

	def hasSequenceDescent(self):
		if self.sequence != 'N':
			return True
		elif:
			any(X.hasSequenceDescent() for X in self.children)
	
	def isSequenceJunction(self):
		return sum(X.hasSequenceDescent() for X in self.children) > 1

	def isSequenceAmbiguous(self):
		return self.sequence == 'N' and self.isSequenceJunction()

	def leftSide(self):
		return Side(self, True)

	def rightSide(self):
		return Side(self, False)

	def isBreakendAmbiguous(self):
		return self.leftSide().isAmbiguous() or self.rightSide().isAmbiguous()

	def isAmbiguous(self):
		return self.isTreeAmbiguous() or self.isSequenceAmbiguous() or self.isBreakendAmbiguous()

	def branchSequence(self):
		pass
				
	def sequenceComplexity(self):
		return sum(X.branchSequence() != self.sequence for X in self.children)

	def sequenceAncestor2(self):
		if self.sequence != 'N' or self.parent is None:
			return self
		else:
			return sequenceAncestor2(self.parent)

	def sequenceAncestor(self):
		return sequenceAncestor2(self.parent)
