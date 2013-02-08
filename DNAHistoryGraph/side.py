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
			return self.segment.parent.right
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
		return len(self.children()) >= 2 and sum(X._hasAttachedDescent() for X in self.children()) >= 2

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

	def _liftedBonds2(self,):
		if self.bond is None:
			return self.liftedBonds()
		else:
			return set([ self ])
	
	def liftedBonds(self):
		""" Returns set of labeled segments whose lifting ancestor is self"""
		if len(self.children()) == 0:
			return set()
		return set(reduce(lambda x, y : x | y, [x._liftedBonds2() for x in self.children() ]))
	
	def nonTrivialLiftedBonds(self):
		""" Return list of non trivial lifted edge partners """
		def fn(x):
			#Returns True if no unattached junction on path to ancestor, else false.
			if x.parent() == x.ancestor():
				return False
			if x.parent().bond == None and x.parent().isJunction():
				return True
			return fn(x.parent())
		return set([ x for x in self.liftedBonds() if x.bond.ancestor() != self.bond 
				or fn(x) or fn(x.bond) ])

	##############################
	## Ambiguity
	##############################

	def rearrangementAmbiguity(self):
		# Note: self looping lifted edges are already reported twice so the formula is correct
		if self.bond == None and self.parent() != None:
			return 0
		return max(0, len(self.nonTrivialLiftedBonds()) - 1)
	
	##############################
	## Modules
	##############################
	
	def isModuleMaterial(self):
		return self.bond != None or (self.parent() == None and len(self.liftedBonds()) > 0)

	##############################
	## Output
	##############################

	def dot(self):
		def fn(side):
			if side.left:
				return "normal"
			return "inv"
		l = []
		if self.bond is not None: 
			if self <= self.bond:
				l.append("%i -> %i [color=red, dir=both, arrowtail=%s, arrowhead=%s]" % (id(self.segment), id(self.bond.segment), fn(self), fn(self.bond)))
		if self.isModuleMaterial():
			for descendant in self.nonTrivialLiftedBonds():
				linkedAncestor = descendant.bond.ancestor()
				if self < linkedAncestor:
					l.append("%i -> %i [color=magenta, dir=both, arrowtail=%s, arrowhead=%s]" % (id(self.segment), id(linkedAncestor.segment), fn(self), fn(linkedAncestor)))
		return "\n".join(l)

	##############################
	## Validation
	##############################
	def validate(self):
		assert self.bond is None or self.bond.bond is self
		assert self.opposite is not None and self.opposite.opposite is self
		return True
