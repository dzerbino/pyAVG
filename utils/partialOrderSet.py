#!/sur/bin/env python

from exceptions import RuntimeError

"""Definition of Partial Order Set"""

class PartialOrderSet(set):
	"""Partially ordered set (poset) with online global ordering of its elements"""

	###################################
	## Basics
	###################################
	def __init__(self, iter=[]):
		""" Creates an unconstrained Poset with the elements of iter (if provided) """
		super(PartialOrderSet, self).__init__(iter)
		self.roots = set()
		self.parents = dict()
		self.children = dict()
		self.depth = dict()
		if iter is not None:
			map(lambda X: self.add(X), iter)
	
	def __copy__(self):
		new = PartialOrderSet(self)
		new.roots = copy.copy(self.roots)
		new.parents = copy.copy(self.parents)
		new.children = copy.copy(self.children)
		self.depth = copy.copy(self.depth)
		return new

	###################################
	## Adding stuff
	###################################
	def add(self, element):
		"""Adds an unconstrained element to the set"""
		super(PartialOrderSet, self).add(element)
		self.roots.add(element)
		self.parents[element] = set()
		self.children[element] = set()
		self.depth[element] = len(self)

	def _addEdge(self, ancestral, derived):
		"""Adds an order constraint between two elements (low level)"""
		if ancestral is None or derived is None:
			return
		else:
			assert ancestral in self
			assert derived in self
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

	def remove(self, elem):
		"""Removes element from a poset and all the incident constraints"""
		if elem in self.roots:
			self.roots.remove(elem)

		if elem in self.parents:
			for parent in self.parents[elem]:
				self.children[parent].remove(elem)
			del self.parents[elem]

		if elem in self.children:
			for child in self.children[elem]:
				self.parents[child].remove(elem)
			del self.children[elem]

		depth = self.depth[elem]
		del self.depth[elem]
		self.depths = dict((X, self._correctedDepth(X, depth)) for X in self.depth)
		super(PartialOrderSet, self).remove(elem)
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

	def testConstraint(self, ancestral, derived):
		"""
		Tests whether adding an ordering constraint between two elements in the set could be done without creating a cycle
		"""
		return self.depth[derived] >= self.depth[ancestral] or self._dfsForward(derived, self.depth[ancestral]) is not None

	def addConstraint(self, ancestral, derived):
		"""
		Adds an ordering constraint between two elements in the set, updating the ordering if necessary. 
		Refuses the addition and raises RuntimeError if a contradiction would be created by the addition.
		"""
		assert ancestral in self
		assert derived in self

		lower = self.depth[derived]
		upper = self.depth[ancestral]
		if lower < upper:
			# Apparent contradiction
			RForward = self._dfsForward(derived, upper)
			if RForward is None:
				# Oops, just created a self loop
				raise RuntimeError
				assert False
			else:
				RBackward = self._dfsBackward(ancestral, lower)
				self._reassign(RForward, RBackward)
				self._addEdge(ancestral, derived)
		else:
			# No need to change anything
			self._addEdge(ancestral, derived)

	################################################
	## Removal of constraint
	################################################
	def removeConstraint(self, parent, child):
		"""Removes an ordering constraint between two elements in the set"""
		self.parents[child].remove(parent)
		if len(self.parents[child]) == 0:
			self.roots.add(child)

	################################################
	## Validate
	################################################
	def _validateChildren(self):
		assert all(X in self for X in self.children)
		for X in self.children:
			for Y in self.children[X]:
				Y in self
				X in self.parents[Y]
				self.depth[Y] > self.depth[X]

	def _validateParents(self):
		assert all(X in self for X in self.parents)
		for X in self.parents:
			for Y in self.parents[X]:
				Y in self
				X in self.children[Y]
				self.depth[Y] < self.depth[X]

	def _validateRoots(self):
		assert all(X in self for X in self.roots)
		assert all(len(self.parents[X]) == 0 for X in self.roots)

	def _validateDepth(self):
		assert all(X in self for X in self.depth)

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
	def dot2(self, elem):
		return "\n".join(["%i -> %i" % (id(elem), id(child)) for child in self.children[elem]])

	def dot(self):
		return "\n".join(["digraph G {"] + [self.dot2(X) for X in self.children] + ["}"]) 
###########################################
## Unit test
###########################################
def test_main():
	pos = PartialOrderSet()
	pos.add(2)
	pos.add(3)
	pos.add(1)
	pos.addConstraint(1,2)
	pos.addConstraint(2,3)
	pos.addConstraint(1,3)
	try:
		pos.addConstraint(3,1)
	except RuntimeError:
		assert pos.validate()
		return
	assert False
	
if __name__ == "__main__":
	test_main()
