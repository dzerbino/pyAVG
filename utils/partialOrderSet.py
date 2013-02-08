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

		for parent in self.parents[elem]:
			self.children[parent].remove(elem)
		del self.parents[elem]

		for child in self.children[elem]:
			assert child in self.parents
			assert elem is not child
			assert elem in self.parents[child]
			self.parents[child].remove(elem)
		del self.children[elem]

		depth = self.depth[elem]
		del self.depth[elem]
		self.depth = dict((X, self._correctedDepth(X, depth)) for X in self.depth)
		super(PartialOrderSet, self).remove(elem)
	
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
	
	def _sortF(self, vals):
		# The use of set() to remove duplicates becomes crucial when interleaving the two lists in _reassign
		return sorted(list(set(vals)), key = lambda X: self.depth[X])

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
		if ancestral is derived:
			raise RuntimeError

		lower = self.depth[derived]
		upper = self.depth[ancestral]
		if lower < upper:
			# Apparent contradiction
			RForward = self._dfsForward(derived, upper)
			if RForward is None:
				# Oops, just created a self loop
				raise RuntimeError
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
		self.children[parent].remove(child)
		if len(self.parents[child]) == 0:
			self.roots.add(child)
			
	################################################
	## Checking ancestry
	################################################
	def _ancestors(self, elem):
		# Yes, I know, recursion would be nicer, but Python is shite with deep recursions
		todo = [elem]
		ancestors = list()
		while len(todo) > 0:
			elem = todo.pop()
			if elem not in ancestors:
				ancestors.extend(self.parents[elem])
				todo.extend(self.parents[elem])
		return ancestors

	def _isAncestor(self, parent, child):
		if self.depth[parent] >= self.depth[child]:
			return False
		else:
			return parent in self._ancestors(child)

	def compare(self, elemA, elemB):
		""" Return -1 if elemA is ancestor of elemB, return 1 is elemA is descendant of elemB else return 0 """
		if self._isAncestor(elemA, elemB):
			return -1
		if self._isAncestor(elemB, elemA):
			return 1
		else:
			return 0

	################################################
	## Checking ancestry
	################################################
	def _isAncestor(self, parent, child):
		return not self.testConstraint(child, parent)

	def compare(self, elemA, elemB):
		""" Return -1 if elemA is ancestor of elemB, return 1 is elemA is descendant of elemB else return 0 """
		if self._isAncestor(elemA, elemB):
			return -1
		if self._isAncestor(elemB, elemA):
			return 1
		else:
			return 0

	################################################
	## Validate
	################################################
	def _validateChildren(self):
		assert all(X in self for X in self.children)
		for X in self.children:
			for Y in self.children[X]:
				assert Y in self
				assert X in self.parents[Y]
				assert self.depth[Y] > self.depth[X], self.dot()
		return True

	def _validateParents(self):
		assert all(X in self for X in self.parents)
		for X in self.parents:
			for Y in self.parents[X]:
				assert Y in self
				assert X in self.children[Y]
				assert self.depth[Y] < self.depth[X]
		return True

	def _validateRoots(self):
		assert all(X in self for X in self.roots)
		assert all(len(self.parents[X]) == 0 for X in self.roots)
		return True

	def _validateDepth(self):
		assert all(X in self for X in self.depth)
		assert len(self) == len(self.depth)
		assert len(set(self.depth.values())) == len(self.depth.values()), str(sorted(self.depth.values())) + "\n" + self.dot()
		return True

	def validate(self):
		"""
		Validation
		"""
		assert self._validateChildren()
		assert self._validateParents()
		assert self._validateDepth()
		assert self._validateRoots()
		return True

	###########################################
	## Unit test
	###########################################
	def dot2(self, elem):
		return "\n".join(["%i -> %i" % (id(elem), id(child)) for child in self.children[elem]])

	def depthDot(self, elem):
		return '\n'.join(['%i [label="%i (%i)"]' % (id(elem), id(elem), self.depth[elem])])

	def dot(self):
		return "\n".join(["digraph G {"] + [self.depthDot(X) for X in self.depth] + [self.dot2(X) for X in self.children] + ["}"]) 

###########################################
## Unit test
###########################################
def test_main():
	pos = PartialOrderSet()
	pos.add(2)
	pos.add(3)
	pos.add(1)
	pos.add(4)
	pos.add(5)
	pos.addConstraint(1,2)
	pos.addConstraint(2,3)
	pos.addConstraint(1,3)
	pos.addConstraint(1,4)
	try:
		pos.addConstraint(3,1)
	except RuntimeError:
		assert pos.validate()
		assert pos.compare(1,2) == -1
		assert pos.compare(1,4) == -1
		assert pos.compare(3,2) == 1
		assert pos.compare(1,1) == 0
		assert pos.compare(2,5) == 0
		return
	assert False
	
if __name__ == "__main__":
	test_main()
