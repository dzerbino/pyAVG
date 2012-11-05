#!/usr/bin/env python

import side
import thread
from traversal import Traversal
from label import Label

class Segment(object):
	""" DNA history segment """

	##########################
	## Basics
	##########################
	def __init__(self, sequence=None, parent = None, children = []):
		if sequence is not None:
			self.label = Label(sequence)
		else:
			self.label = None
		self.parent = parent
		self.children = set(children)
		self.right = None
		self.left = side.Side(self, True)
		self.right = side.Side(self, False)

	def __cmp__(self, other):
		return cmp(id(self), id(other))

	def __copy__(self):
		return Segment(self.label)

	def __str__(self):
		return str(self.label)

	def __hash__(self):
		return id(self)

	def sides(self):
		""" Returns list of segment sides """
		return [self.left, self.right]

	def createBranch(self, other):
		""" Creates branch between segments """
		self.children.add(other)
		other.parent = self

	def deleteBranch(self, other):
		""" Removes branch between segments """
		self.children.remove(other)
		other.parent = None

	def getSide(self, left):
		if left:
			return self.left
		else:
			return self.right

	##########################
	## Lifted labels
	##########################

	def disconnect(self):
		"""Destroy pointers to this segment """
		self.left.deleteBond()
		self.right.deleteBond()
		for child in self.children:
			child.parent = self.parent
		if self.parent is not None:
			self.parent.children |= self.children
			self.parent.children.remove(self)

	def _ancestor2(self):
		if self.label is not None or self.parent is None:
			return self
		else:
			return self.parent._ancestor2()

	def ancestor(self):
		if self.parent is None:
			return self
		else:
			return self.parent._ancestor2()

	def _liftedLabels2(self):
		if self.label is None:
			liftingLabels = sum([X._liftedLabels2() for X in self.children], [])
			if len(liftingLabels) < 2:
				return liftingLabels
			else:
				return [(X[0], True) for X in liftingLabels]
		else:
			return [(self.label, self.label != self.ancestor().label)]
	
	def liftedLabels(self):
		""" Returns tuple of lifted labels (X, Y), where X is lifted label, and Y determines whether X is non-trivial """
		return sum([X._liftedLabels2() for X in self.children], [])
	
	def isJunction(self):
		return sum(X._liftedLabels2() > 0 for X in self.children) > 1

	def nonTrivialLiftedLabels(self):
		if self.label is not None:
			return [X[0] for X in self.liftedLabels() if X[1]]
		else:
			return [X[0] for X in self.liftedLabels()]

	##########################
	## Ambiguity
	##########################
	def coalescenceAmbiguity(self):
		return max(0, len(self.children) - 2)

	def substitutionAmbiguity(self):
		return max(0, len(self.nonTrivialLiftedLabels()) - 1)

	def rearrangementAmbiguity(self):
		return self.left.rearrangementAmbiguity() + self.right.rearrangementAmbiguity()

	##########################
	## Cost
	##########################
	def substitutionCost(self, lowerBound = True):
		if lowerBound:
			return len(set(self.nonTrivialLiftedLabels()))
		else:
			return len(self.nonTrivialLiftedLabels())

	##########################
	## Threads
	##########################
	def thread(self):
		return thread.Thread([Traversal(self, True)])

	def threads(self, data):
		threads, segmentThreads = data
		if self in segmentThreads:
			return data
		else:
			thread = self.thread()
			for traversal in thread:
				segmentThreads[traversal.segment] = thread
			threads.add(thread)

			return threads, segmentThreads

	##########################
	## Output
	##########################
	def dot(self):
		label = str(self.label)	
		if self.substitutionAmbiguity() > 0:
			label += 'S'
		if self.coalescenceAmbiguity() > 0:
			label += 'C'
		if self.left.rearrangementAmbiguity() > 0:
			label += 'L'
		if self.right.rearrangementAmbiguity() > 0:
			label += 'R'
		lines = ['%i [label="%s"]' % (id(self), label)]
		if self.parent is not None:
			if self.parent.label == self.label:
				lines.append('%i -> %i [color=green]' % (id(self.parent), id(self)))
			else:
				lines.append('%i -> %i [color=blue]' % (id(self.parent), id(self)))
		lines.append(self.left.dot())
		if self.left.bond is not self.right:
			lines.append(self.right.dot())
		return "\n".join(lines)	
		
	
	##########################
	## Validation
	##########################
	def validate(self):
		assert self.parent is None or self in self.parent.children
		assert all(self is child.parent for child in self.children)
		assert self.left.validate()
		assert self.right.validate()
		assert self.substitutionCost() <= self.substitutionCost(False)
		return True
