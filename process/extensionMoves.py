#!/usr/bin/env python

import random
from pyAVG.DNAHistoryGraph.label import Label

class ExtensionMove(object):
	""" Wrapper for possible graph transformation """
	def __init__(self, function, args):
		self.args = args
		self.function = function

def getMajority(list):
	counters = dict()
	for x in list:
		if x in counters:
			counters[x] += 1
		else:
			counters[x] = 1
	max = -1
	winner = None
	for x in counters:
		if counters[x] > max:
			max = counters[x]
			winner = x
	return winner

###############################################
## Case 1
## Adding necessary labels
###############################################

def listCase1(graph):
	return [ExtensionMove(applyCase1, (segment, graph)) for segment in filter(lambda X: X.substitutionAmbiguity() > 0, graph.segments)]

def applyCase1(args):
	segment, graph = args
	if segment.label is None:
		# In an ideal world Fitch parsimony would be a nice touch...
		segment.label = Label(random.choice(['A','T'])) 
	else:
		children = list(segment.children)
		children.pop(random.randrange(len(children)))
		if len(children) == 1 and children[0].label is None and random.random() < 0.5:
			bridge = children[0]
		else:
			bridge = graph.newSegment()
			for child in children:
				graph.deleteBranch(segment, child)
				graph.createBranch(bridge, child) 
			graph.createBranch(segment, bridge)
		bridge.label = Label(segment.label)
	graph.validate()

###############################################
## Case 2
## Adding necessary bonds
###############################################

def listCase2(graph):
	return [ExtensionMove(applyCase2, (side, graph)) for side in filter(lambda X: X.rearrangementAmbiguity() > 0, graph.sides())]

def applyCase2(args):
	side, graph = args
	if side.bond is None:
		# In an ideal world Fitch parsimony would be a nice touch...
		ancestor = side.ancestor()
		if ancestor is not None and ancestor.bond is not None:
			# Pulling down
			# Note: not ping pong
			segment = graph.newSegment()
			graph.createBranch(ancestor.bond.segment, segment)
			graph.validate()
			assert graph.eventGraph.testConstraint(graph.sideThread(ancestor), graph.sideThread(side))
			assert graph.eventGraph.testConstraint(graph.sideThread(ancestor.bond), graph.sideThread(side))
			graph.createBond(side, segment.getSide(ancestor.bond.left))
			graph.validate()
		else:	
			# Pulling up
			liftedEdges = [X[0] for X in side.liftedPartners()]
			target = getMajority(liftedEdges)

			if graph.eventGraph.testConstraint(graph.sideThread(target), graph.sideThread(side)):
				if graph.eventGraph.testConstraint(graph.sideThread(side), graph.sideThread(target)):
					# side must be connected to that very target, not a child
					if target.bond is not None:	
						# Stubbing
						segment = graph.newSegment()
						target = segment.left
						graph.createBond(side, target)
						graph.validate()
					else:
						graph.createBond(side, target)
						graph.validate()
				else:
					# Pulling down from target
					segment = graph.newSegment()
					graph.createBranch(target.segment, segment)
					graph.createBond(side, segment.getSide(target.left))
					graph.validate()
			else:
				# Stubbing
				segment = graph.newSegment()
				graph.createBond(side, segment.left)
				graph.validate()
			
	else:
		# Bridging
		children = list(side.children())
		children.pop(random.randrange(len(children)))
		if len(children) == 1 and children[0].bond is None and random.random() < 0.5:
			bridge = children[0]
		else:
			bridgeS = graph.newSegment()
			for child in children:
				graph.deleteBranch(side.segment, child.segment)
				graph.createBranch(bridgeS, child.segment) 
				graph.validate()
			graph.createBranch(side.segment, bridgeS)
			graph.validate()
			bridge = bridgeS.getSide(side.left)
		# Note: created segment is hanging but the ancestor of its partner, i.e. side, is not hanging, ergo no ping pong
		segment = graph.newSegment()
		graph.createBranch(side.bond.segment, segment)
		graph.validate()
		graph.createBond(bridge, segment.getSide(side.bond.left))
		graph.validate()
		
###############################################
## Case 3
## Adding necessary coalescences
###############################################

def listCase3(graph):
	return [ExtensionMove(applyCase3, (segment, graph)) for segment in filter(lambda X: X.coalescenceAmbiguity() > 0, graph.segments)]

def applyCase3(args):
	segment, graph = args 
	children = list(segment.children)
	children.pop(random.randrange(len(children)))	
	bridge = graph.newSegment()
	for child in children:
		graph.deleteBranch(segment, child)
		graph.createBranch(bridge, child) 
	graph.createBranch(segment, bridge)
	graph.validate()

detectors = [listCase1, listCase2, listCase3]
