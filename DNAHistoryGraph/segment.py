#!/usr/bin/env python

import side
from thread import Thread

class Segment(object):
	""" DNA history segment """

	##########################
	## Basics
	##########################
	def __init__(self, sequence=None, parent = None, children = []):
		self.sequence = str(sequence)[:1]
		assert parent is None or isinstance(parent, Segment)
		self.parent = parent
		self.children = set(children)
		self.right = None
		self.left = side.Side(self, True)
		self.right = side.Side(self, False)

	def __cmp__(self, other):
		return cmp(id(self), id(other))

	##########################
	## Lifted labels
	##########################

	def _ancestor2(self):
		if self.sequence != None or self.parent is None:
			return self
		else:
			return self.parent._ancestor2()

	def ancestor(self):
		if self.parent is None:
			return self
		else:
			return self.parent._ancestor2()

	def _liftedLabels2(self):
		if self.sequence is None:
			liftingLabels = sum([X._liftedLabels2() for X in self.children], [])
			if len(liftingLabels) < 2:
				return liftingLabels
			else:
				return [(X[0], True) for X in liftingLabels]
		else:
			return [(self.sequence, self.sequence is not self.ancestor().sequence)]
	
	def liftedLabels(self):
		""" Returns tuple of lifted labels (X, Y), where X is lifted label, and Y determines whether X is non-trivial """
		return sum([X._liftedLabels2() for X in self.children], [])

	def nonTrivialLiftedLabels(self):
		return [X[0] for X in self.liftedLabels() if X[1]]

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
	def expandThread(self, thread):
		if self in thread:
			return thread
		else:
			thread.append(self)
			return self.left.expandThread(self.right.expandThread(thread))

	def threads(self, data):
		threads, segmentThreads = data
		if self in segmentThreads:
			return data
		else:
			thread = self.expandThread(Thread())
			for segment in thread:
				segmentThreads[segment] = thread
			threads.add(thread)

			return threads, segmentThreads

	##########################
	## Modules
	##########################
	def modules(self, data):
		print 'SEG'
		return reduce(lambda X, Y: Y.modules(X), [self.left, self.right], data)

	##########################
	## Output
	##########################
	def dot(self):
		node = " ".join([str(id(self)), "[ label=", self.sequence, "]"])
		if self.parent is None:
			parent = " ".join([str(id(self.parent)), "->", str(id(self))])
		else:
			parent = ""
		if self.left.bond is not None and self < self.left.bond.segment:
			left = " ".join([str(id(self)), "->", str(id(self.left.bond.segment))])
		else:
			left = ""
		if self.right.bond is not None and self < self.right.bond.segment:
			right = " ".join([str(id(self)), "->", str(id(self.right.bond.segment))])
		else:
			right = ""
		return "\n".join([node, parent, left, right])	
		
	
	##########################
	## Validation
	##########################
	def validate(self):
		assert self.parent is None or self in self.parent.children
		assert all(self is child.parent for child in self.children)
		assert self.left.validate()
		assert self.right.validate()
		return True
