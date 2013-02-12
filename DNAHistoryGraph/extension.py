#!/usr/bin/env python

import copy
from graph import DNAHistoryGraph
from label import Label

class GraphExtension(DNAHistoryGraph):
	def __init__(self, reduction):
		duplicates = dict(map(lambda X: (X, copy.copy(X)), reduction.segments))
		for segment in reduction.segments:
			if segment.left.bond is not None and duplicates[segment].left.bond is None:
				duplicates[segment].left.createBond(duplicates[segment.left.bond.segment].getSide(segment.left.bond.left))
			if segment.right.bond is not None and duplicates[segment].right.bond is None:
				duplicates[segment].right.createBond(duplicates[segment.right.bond.segment].getSide(segment.right.bond.left))
			if segment.parent is not None and duplicates[segment].parent is None:
				duplicates[segment.parent].createBranch(duplicates[segment])
		super(GraphExtension, self).__init__(duplicates)
		self.irreducibleSegments = frozenset(self.segments) 
		self.irreducibleSides = frozenset(filter(lambda X: X.bond is not None, self.sides()))
		self.irreducibleLabels = frozenset(segment.label for segment in self.segments if segment.label is not None)

	def isGBounded(self):
		return hasNoGReducibleLabels(self) and hasNoGReducibleSegments(self) and hasNoGReducibleBonds(self) and hasNoGReduciblePingPongs(self)

	def makeGBounded(self):
		while self.removeGReducibleElements():
			pass

	def removeGReducibleElements(self):
		A = removeGReducibleLabels(self) 
		B = removeGReducibleSegments(self) 
		C = removeGReducibleBonds(self)
		D = removeGReduciblePingPongs(self)
		return A > 0 or B > 0 or C > 0 or D > 0
	
	##################################
	#Remove non-minimal elements
	##################################
	
	def removeNonMinimalElements(self):	
		for segment in list(self.segments):
			if segment.label != None and (not segment.isLabelJunction()) and (not segment.isLabelBridge()):
				segment.deleteLabel()
			
			for side in (segment.left, segment.right): #Get rid of useless bonds
				if side.bond != None and (not side.isJunction()) and (not side.bond.isJunction()) and (not side.isBridge()):
					self.deleteBond(side)
		
		for segment in list(self.segments): #Get rid of useless nodes
			if segment.label == None and segment.left.bond == None and segment.right.bond == None:
				segment.disconnect()
				self.segments.remove(segment) 
		
		# Recompute the event graph from scratch
		self.eventGraph, self.segmentThreads = self.threads()
		self.timeEventGraph()

	##################################
	## Output
	##################################
	def irreducibleDot(self):
		return ['node [style=filled]'] + [str(id(X)) for X in self.irreducibleSegments] + ['node [style=none]']

	def dot(self):
		return "\n".join(["digraph G {"] + self.irreducibleDot() + [X.dot() for X in self.eventGraph] + ["}"])

def isNonTrivialLift(segment):
	ancestor = segment.ancestor()
	if segment.label != ancestor.label:
		return True
	parent = segment.parent
	while parent is not None and parent is not ancestor:
		if parent.isJunction():
			return True
		parent = parent.parent
	return False

def hasGReducibleLabel(segment, graph):
	if segment.label is None or segment.label in graph.irreducibleLabels:
		return False

	L = segment.liftedLabels()
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
	assert not any(hasGReducibleLabel(segment, graph) for segment in graph.segments)
	return True

def removeGReducibleLabel(segment, graph):
	if hasGReducibleLabel(segment, graph):
		segment.label = None
		print 'G-Reducible label REMOVED'
		return 1
	else:
		return 0

def removeGReducibleLabels(graph):
	return sum(map(lambda X: removeGReducibleLabel(X, graph), graph.segments))

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
		print 'LEAFL', id(side.segment)
		return True
	elif (len(liftedA) > 1 or len(liftedB) > 1):
		# Junction or complex
		return False
	elif len(liftedA) == 1 and len(liftedB) == 1 and not liftedA[0][1] and not liftedB[0][1]:
		# Redundant
		print 'REDUNDANT', id(side.segment)
		return True
	elif len(liftedA) == 1 and len(liftedB) == 1 and isNonTrivialLiftingBond(side):
		# Complicating
		print 'COMPLICATING', id(side.segment)
		return True
	else:
		# Bridge
		return False

def hasNoGReducibleBonds(graph):
	assert not any(isGReducibleBond(side, graph) for side in graph.sides())
	return True

def removeGReducibleBond(side, graph):
	if isGReducibleBond(side, graph):
		graph.deleteBond(side)
		print 'G-reducible bond REMOVED'
		return 1
	else:
		return 0

def removeGReducibleBonds(graph):
	return sum(map(lambda X: removeGReducibleBond(X, graph), list(graph.sides())))

def isGReducibleSegment(segment, graph):
	if segment in graph.irreducibleSegments:
		return False
	elif segment.parent is None and len(segment.children) == 0 and segment.left.bond is None and segment.right.bond is None:
		# Isolated segment
		print 'ISOLATED'
		return True
	elif segment.label is None and segment.left.bond is None and segment.right.bond is None and len(segment.children) <= 1:
		# Free tailed branch with at most one child or free headed branch
		print 'FREE TAILED BRANCH W/ > 1 CHILD OR FREE HEADED BRANCH'
		return True
	else:
		# None of the above
		return False
	
def hasNoGReducibleSegments(graph):
	assert not any(isGReducibleSegment(segment, graph) for segment in graph.segments)
	return True

def removeGReducibleSegment(segment, graph):
	if isGReducibleSegment(segment, graph):
		for child in list(segment.children):
			graph.deleteBranch(segment, child)
			if segment.parent is not None:
				graph.createBranch(segment.parent, child)
		if segment.parent is not None:
			segment.parent.children.remove(segment)
		graph.eventGraph.remove(graph.segmentThreads[segment])
		del graph.segmentThreads[segment]
		graph.segments.remove(segment)
		print 'G-Reducible segment REMOVED', id(segment)
		return 1
	else:
		return 0

def removeGReducibleSegments(graph):
	return sum(map(lambda X: removeGReducibleSegment(X, graph), graph.segments))

def isHanging(side):
	return len(side.liftedPartners()) == 0

def isPing(side, graph):
	return side.bond is not None and side not in graph.irreducibleSides and isHanging(side) and side.ancestor().bond is not None and isHanging(side.ancestor().bond)

def hasNoGReduciblePingPongs(graph):
	assert all(not isPing(X, graph) for X in graph.sides())
	return True

def removePingPong(side, graph):
	if isPing(side, graph):
		print 'Ping', id(side.segment)
		exPartner = side.bond
		graph.deleteBond(side)
		ancestor = exPartner.ancestor()
		if exPartner.parent is not None and ancestor.bond is not None:
			# Pull down correction
			new = graph.newSegment()
			graph.createBranch(ancestor.bond.segment, new)
			graph.createBond(exPartner, new.getSide(ancestor.bond.left))
			print 'Ping Pong corrected w/ pull down', id(side.segment)
		else:
			# Stub
			graph.createBond(exPartner, graph.newSegment().left)
			print 'Ping Pong corrected w/ pull stub', id(side.segment)
		return 1
	return 0

def removeGReduciblePingPongs(graph):
	return sum(map(lambda X: removePingPong(X, graph), graph.sides()))
