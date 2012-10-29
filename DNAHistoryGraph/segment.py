#!/usr/bin/env python

from side import Side

class Segment(object):
	""" DNA history segment """

	##########################
	## Basics
	##########################
	def __init__(self, sequence=None, parent = None, children = []):
		self.sequence = str(sequence)[:1]
		assert isinstance(parent, Segment)
		self.parent = parent
		self.children = set(children)
		self.left = Side(self, True)
		self.right = Side(self, False)

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
		if self in thread.segments:
			return thread
		else:
			thread.segments.add(self)
			return self.left.expandThread(self.right.expandThread(thread))

	def threads(self, data):
		threads, segmentThreads = data:
		if self in segmentThreads:
			return data
		else:
			thread = self._expandThread(Thread())
			for segment in thread:
				segmentThreads[segment] = thread
			threads.add(thread)

			return threads, segmentThreads

	##########################
	## Modules
	##########################
	def modules(self, data):
		return reduce(lambda X, Y: Y.modules(X), [self.left, self.right], data)

	##########################
	## Validation
	##########################
	def validate(self):
		assert self.parent is None or self in self.parent.children
		assert all(self is child.parent for child in self.children)
		assert self.left.validate()
		assert self.right.validate()
