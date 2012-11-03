#!/usr/bin/env python

import copy
from graph import DNAHistoryGraph

class GraphExtension(DNAHistoryGraph):
	def __init__(self, reduction):
		duplicates = dict(map(lambda X: (X, copy.copy(X)), reduction.segments))
		for segment in reduction.segments:
			if segment.left.bond is not None:
				if segment.left.bond.left:
					duplicates[segment].left.createBond(duplicates[segment.left.bond.segment].left)
				else:
					duplicates[segment].left.createBond(duplicates[segment.left.bond.segment].right)
			if segment.right.bond is not None:
				if segment.right.bond.left:
					duplicates[segment].right.createBond(duplicates[segment.right.bond.segment].left)
				else:
					duplicates[segment].right.createBond(duplicates[segment.right.bond.segment].right)

			if segment.parent is not None:
				duplicates[segment].parent = duplicates[segment.parent]
			duplicates[segment].children = set(duplicates[X] for X in segment.children)
		super(GraphExtension, self).__init__(duplicates)
		self.irreducibleSegments = list(self.segments) 
		self.irreducibleSides = filter(lambda X: X.bond is not None, self.sides()) 	
		self.irreducibleLabels = [segment.label for segment in self.segments if segment.label is not None]

	def isGBounded(self):
		return hasNoGReducibleLabels(self) and hasNoGReducibleSegments(self) and hasNoGReducibleBonds(self) and hasNoGReduciblePingPongs(self)

def isNonTrivialLift(segment):
	ancestor = segment.ancestor()
	if segment.label != ancestor.label:
		return True
	parent = segment.parent
	while parent is not ancestor:
		if parent.isJunction():
			return True
		parent = parent.parent
	return False

def hasGReducibleLabel(segment, graph):
	if segment.label is None or segment.label in graph.irreducibleLabels:
		return False

	L = segment.liftedEdges()
	if len(L) == 0:
		# Leaf
		return True
	elif len(L) > 1:
		# junction or complex
		return False
	else:
		# Only one lifted edge
		if L[0][0] == segment.label:
			# Redundant label
			return True
		else:
			if isNonTrivialLift(segment):
				# Complicating label
				return True
			else:
				# Bridge
				return False

def hasNoGReducibleLabels(graph):
	return not any(hasGReducibleLabel(segment, graph) for segment in graph.segments)

def isNonTrivialLiftingBond(side):
	ancestorA = side.ancestor()
	ancestorB = side.bond.ancestor()	
	if ancestorB is not ancestorA.bond:
		return True

	if ancestorA is not side:
		parent = side.parent()
		while parent is not ancestorA:
			if parent.isJunction():
				return True
			parent = parent.parent()

	if ancestorB is not side.bond:
		parent = side.bond.parent()
		while parent is not ancestorB:
			if parent.isJunction():
				return True
			parent = parent.parent()

	return False

def isGReducibleBond(side, graph):
	if side.bond is None or side in graph.irreducibleSides:
		return False

	# To avoid double analysis
	if side.bond < side:
		return False

	liftedA = side.liftedPartners()		
	liftedB = side.bond.liftedPartners()

	if len(liftedA) == 0 and len(liftedB) == 0:
		# Leaf
		return True
	elif (len(liftedA) > 1 or len(liftedB) > 1):
		# Junction or complex
		return False
	elif len(liftedA) == 1 and len(liftedB) == 1 and not liftedA[0][1] and not liftedB[0][1]:
		# Redundant
		return True
	elif isNonTrivialLiftingBond(side):
		# Complicating
		return True
	else:
		# Bridge
		return False

def hasNoGReducibleBonds(graph):
	return not any(isGReducibleBond(side, graph) for side in graph.sides())

def isGReducibleSegment(segment, graph):
	if segment in graph.irreducibleSegments:
		return False
	elif segment.parent is None and len(segment.children) == 0 and segment.left.bond is None and segment.right.bond is None:
		# Isolated segment
		return True
	elif segment.label is None and segment.left.bond is None and segment.right.bond is None and segment.parent is None and len(segment.children) == 1:
		# Free headed branch
		return True
	elif segment.label is None and segment.left.bond is None and segment.right.bond is None and len(segment.children) == 0:
		# Free tailed branch with at most one child
		return True
	else:
		# None of the above
		return False
	
def hasNoGReducibleSegments(graph):
	return not any(isGReducibleSegment(segment, graph) for segment in graph.segments)

def isHanging(side):
	return len(side.liftedPartners()) == 0

def isPing(side):
	return isHanging(side) and isHanging(side.ancestor().bond)

def hasNoGReduciblePingPongs(graph):
	return all(lambda X: not isPing(X) for X in graph.sides())
