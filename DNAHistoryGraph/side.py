#!/usr/bin/env python

class Side(object):
	""" Segment side in DNA history graph """

	##############################
	## Basics
	##############################
	def __init__(self, segment, left):
		self.segment = segment
		self.left = bool(left)
		self.bond = None
		if self.left:
			self.opposite = self.segment.right
		else:
			self.opposite = self.segment.left
		if self.opposite is not None:
			self.opposite.opposite = self

	def createBond(self, other):
		self.bond = other
		other.bond = self

	def deleteBond(self):
		if self.bond is not None:
			self.bond.bond = None
			self.bond = None

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
			return [X.left for X in self.segment.children]
		else:
			return [X.right for X in self.segment.children]

	##############################
	## Lifted edges
	##############################
	def _ancestor2(self):
		if self.bond is not None or self.parent() is None:
			return self
		else:
			return self.parent()._ancestor2()

	def ancestor(self):
		if self.parent() is None:
			return self
		else:
			return self.parent()._ancestor2()

	def _liftedBonds2(self):
		if self.bond is None:
			liftingBonds = sum([X._liftedBonds2() for X in self.children()], [])
			if len(liftingBonds) < 2:
				return liftingBonds
			else:
				return [(X[0], True) for X in liftingBonds]
		else:
			target = self.bond.ancestor()
			return [(target, target is not self.ancestor().bond)]

	def liftedBonds(self):
		return sum([X._liftedBonds2() for X in self.children()], [])

	def nonTrivialLiftedBonds(self):
		return [X[0] for X in self.liftedBonds() if X[1]]

	##############################
	## Ambiguity
	##############################

	def rearrangementAmbiguity(self):
		return max(0, len(self.nonTrivialLiftedBonds()) - 1)

	##############################
	## Threads
	##############################
	def expandThread(self, thread):
		if self.bond is None:
			return thread
		else:
			return self.bond.segment.expandThread(thread)

	##############################
	## Modules
	##############################
	def _expandModule(self, module):
		module.sides.add(self)
		if self.bond not in module.sides:
			self.bond._expandModule(module)
		for liftedEdge in self.nonTrivialLiftedEdges():
			if liftedEdge not in module.sides():
				module.nonTrivialLiftedEdges.add(frozenset(self, liftedEdge))
				liftedEdge._expandModule(module)

	def modules(self, data):
		""" Returns a list of modules and a set of already visited sides """
		data = modules, visited
		if self in visited:
			return data
		else:
			module = self._expandModule(Module())
			return modules + [module], visited + module.sides

	##############################
	## Validation
	##############################
	def validate(self):
		assert self.bond is None or self.bond.bond is self
		assert self.opposite is not None and self.opposite.opposite is self
		return True
