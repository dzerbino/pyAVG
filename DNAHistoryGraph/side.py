#!/usr/bin/env python

import module
import liftedEdge

def _junctionsOnTheWay(side):
	ancestor = side.ancestor()
	parent = side.parent()
	while parent is not None and parent is not ancestor:
		if parent.isJunction():
			return True
		parent = parent.parent()
	return False

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

	def __cmp__(self, other):
		return cmp(id(self), id(other))

	def __hash__(self):
		return id(self)

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
	def _hasAttachedDescent(self):
		if self.bond is not None:
			return True
		else:
			return any(X._hasAttachedDescent() for X in self.children())

	def isJunction(self):
		return len(self.children()) > 2 and sum(X._hasAttachedDescent() for X in self.children()) > 2

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

	def _liftedPartners2(self):
		""" Recursive element of liftedPartners """
		if self.bond is None:
			childLiftedBonds = [X._liftedPartners2() for X in self.children()]
			liftingPartners = sum(childLiftedBonds, [])
			if sum(len(X) > 0 for X in childLiftedBonds) > 1:
				# Unattached bond junction!!
				return [(X[0], True, X[2]) for X in liftingPartners]
			else:
				# Unattached bond on linear lifting path 
				return liftingPartners
		else:
			target = self.bond.ancestor()
			return [(target, target is not self.ancestor().bond or _junctionsOnTheWay(target), self)]

	def liftedBonds(self):
		""" Returns list of tuples (lifted edge partner of self, is non trivial) """
		if self.parent() is not None or self.bond is None:
			return sum([X._liftedPartners2() for X in self.children()], [])
		else:
			# Special case for root nodes which are their own ancestors
			return sum([X._liftedPartners2() for X in self.children()], [(self.bond.ancestor(), self.bond.segment.parent is not None)])

	def nonTrivialLiftedBonds(self):
		""" Return list of non trivial lifted edge partners """
		if self.bond is not None:
			return [X[0] for X in self.liftedPartners() if X[1]]
		else:
			return [X[0] for X in self.liftedPartners()]

	##############################
	## Ambiguity
	##############################

	def rearrangementAmbiguity(self):
		# Note: self looping lifted edges are already reported twice so the formula is correct
		return max(0, len(self.nonTrivialLiftedBonds()) - 1)

	##############################
	## Modules
	##############################
	def _expandModule(self, module):
		selfLoops = 0
		module.sides.add(self)
		if self.bond is not None and self.bond not in module.sides:
			self.bond._expandModule(module)
		for partner in self.nonTrivialLiftedPartners():
			if partner is self:
				selfLoops += 1
			if partner < self:
				module.nonTrivialLiftedEdges.append(liftedEdge.LiftedEdge((self, partner)))
			if partner not in module.sides:
				partner._expandModule(module)
		for i in range(selfLoops / 2):
			# Note: each self loop is reported twice in the list (one for each incidence)	
			module.nonTrivialLiftedEdges.append(liftedEdge.LiftedEdge((self, self)))

	def modules(self, data):
		""" Returns a list of modules and a set of already visited sides """
		modules, visited = data
		if self in visited:
			return data
		else:
			M = module.Module()
			self._expandModule(M)
			return modules + [M], visited | M.sides

	def isModuleMaterial(self):
		return self.bond is not None or (self.segment.parent is None and len(self.liftedPartners()) > 1)

	##############################
	## Output
	##############################

	def dot(self):
		if self.bond is not None and self.segment <= self.bond.segment:
			return "%i -> %i [color=red, arrowhead=none]" % (id(self.segment), id(self.bond.segment)) 
		else:
			return ""

	##############################
	## Validation
	##############################
	def validate(self):
		assert self.bond is None or self.bond.bond is self
		assert self.opposite is not None and self.opposite.opposite is self
		return True
