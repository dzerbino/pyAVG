#!/sur/bin/env python

"""Definition of Partial Order Set"""

class PartialOrderSet(object):
	"""Partially ordered set (poset) with online global ordering of its elements"""

	###################################
	## Basics
	###################################
	def __init__(self, iter=None):
		""" Creates an unconstrained Poset with the elements of iter (if provided) """
		self.roots = set()
		self.elements = list()
		self.parents = dict()
		self.children = dict()
		self.depth = dict()
		if iter is not None:
			map(lambda X: self.addElement(X), iter)
	
	def __copy__(self):
		new = PartialOrderSet(self.elements)
		new.roots = copy.copy(self.roots)
		new.parents = copy.copy(self.parents)
		new.children = copy.copy(self.children)
		self.depth = copy.copy(self.depth)
		return new

	###################################
	## Adding stuff
	###################################
	def addElement(self, element):
		"""Adds an unconstrained element to the set"""
		self.elements.append(element)
		self.roots.add(element)
		self.parents[element] = set()
		self.children[element] = set()
		self.depth[element] = len(self.elements)

	def _addEdge(self, ancestral, derived):
		"""Adds an order constraint between two elements (low level)"""
		if ancestral is None or derived is None:
			return
		else:
			assert ancestral in self.elements
			assert derived in self.elements
			if derived in self.roots:
				self.roots.remove(derived)
			self.parents[derived].add(ancestral)
			self.children[ancestral].add(derived)

	###################################
	## Removing an element
	###################################
	def _correctedDepth(self, elem, depth):
		oldDepth = self.depth[elem]
		if oldDepth <= depth:
			return oldDepth
		else:
			return oldDepth - 1

	def pop(self, elem):
		"""Removes element from a poset and all the incident constraints"""
		if elem in self.roots:
			self.roots.remove(elem)
		else:
			assert elem in self.parents
			for parent in self.parents[elem]:
				self.children[parent].remove(elem)

		if elem in self.children:
			for child in self.children[elem]:
				self.parents[child].remove(elem)

		depth = self.depth[elem]
		del self.depth[elem]
		self.depths = dict((X, self._correctedDepth(X, depth)) for X in self.depth)
		self.elements.remove(elem)
		return elem

	####################################################
	## Adding a constraint
	####################################################
	def _dfsForward(self, x, upper):
		RForward = [x]
		for y in self.children[x]:
			if self.depth[y] == upper:
				return None
			elif y not in RForward and self.depth[y] < upper:
				recursion = self._dfsForward(y, upper)
				if recursion is None:
					return None
				else:
					RForward += recursion
		return RForward

	def _dfsBackward(self, x, lower):
		RBackward = [x]
		for y in self.parents[x]:
			if self.depth[y] == lower:
				return None
			elif y not in RBackward and self.depth[y] > lower:
				recursion = self._dfsBackward(y, lower)
				if recursion is None:
					return None
				else:
					RBackward += recursion
		return RBackward
	
	def _sortF(self, list):
		return sorted(list, key = lambda X: self.depth[X])

	def _reassign(self, RForward, RBackward):
		Lnodes = self._sortF(RBackward) + self._sortF(RForward)
		Ldepths = sorted([self.depth[X] for X in Lnodes])
		for index in range(len(Lnodes)):
			self.depth[Lnodes[index]] = Ldepths[index]

	def addConstraint(self, ancestral, derived):
		"""Adds an ordering constraint between two elements in the set, updating the ordering if necessary. Refuses the addition and returns False if a contradiction would be created by the addition."""
		if ancestral not in self.elements:
			self.addElement(ancestral)
		if derived not in self.elements:
			self.addElement(derived)

		lower = self.depth[derived]
		upper = self.depth[ancestral]
		if lower < upper:
			# Apparent contradiction
			RForward = self._dfsForward(derived, upper)
			if RForward is None:
				# Oops, just created a self loop
				return False
			else:
				RBackward = self._dfsBackward(ancestral, lower)
				self._reassign(RForward, RBackward)
				self._addEdge(ancestral, derived)
				return True
		else:
			# No need to change anything
			self._addEdge(ancestral, derived)
			return True

	################################################
	## Removal of constraint
	################################################
	def removeConstraint(self, parent, child):
		"""Removes an ordering constraint between two elements in the set"""
		self.parents[child].remove(parent)
		if len(self.children[child]) == 0 and len(self.parents[child]) == 0:
			self.pop(child)
		elif len(self.parents[child]) == 0:
			self.roots.add(child)

		if len(self.children[parent]) == 0 and len(self.parents[parent]) == 0:
			self.pop(parent)

	################################################
	## Validate
	################################################
	def _validateChildren(self):
		assert all(X in self.elements for X in self.children)
		for X in self.children:
			for Y in self.children[X]:
				Y in self.elements
				X in self.parents[Y]
				self.depth[Y] > self.depth[X]

	def _validateParents(self):
		assert all(X in self.elements for X in self.parents)
		for X in self.parents:
			for Y in self.parents[X]:
				Y in self.elements
				X in self.children[Y]
				self.depth[Y] < self.depth[X]

	def _validateRoots(self):
		assert all(X in self.elements for X in self.roots)
		assert all(len(self.parents[X]) == 0 for X in self.roots)

	def _validateDepth(self):
		assert all(X in self.elements for X in self.depth)

	def validate(self):
		"""
		Validation
		"""
		self._validateChildren()
		self._validateParents()
		self._validateDepth()
		self._validateRoots()
		return True

###########################################
## Unit test
###########################################
def test_main():
	pos = PartialOrderSet()
	pos.addElement(2)
	pos.addElement(3)
	pos.addElement(1)
	assert pos.addConstraint(1,2) == True
	assert pos.addConstraint(2,3) == True
	assert pos.addConstraint(1,3) == True
	assert pos.addConstraint(3,1) == False
	assert pos.validate()
	
if __name__ == "__main__":
	main()
