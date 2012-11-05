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
		print 'Adding junction label'
		segment.label = Label(random.choice(['A','T'])) 
	else:
		children = filter(lambda X: X.label is not None or len(X.liftedLabels()) > 0, segment.children)
		children.pop(random.randrange(len(children)))
		if len(children) == 1 and children[0].label is None and random.random() < 0.5:
			print 'Labelling child'
			bridge = children[0]
		else:
			print 'Adding label bridge'
			bridge = graph.newSegment()
			for child in children:
				graph.deleteBranch(segment, child)
				graph.createBranch(bridge, child) 
			graph.createBranch(segment, bridge)
		bridge.label = Label(segment.label)

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
			print 'Pulling down junction bond'
			# Pulling down
			# Note: not ping pong
			segment = graph.newSegment()
			graph.createBranch(ancestor.bond.segment, segment)
			assert graph.eventGraph.testConstraint(graph.sideThread(ancestor), graph.sideThread(side))
			assert graph.eventGraph.testConstraint(graph.sideThread(ancestor.bond), graph.sideThread(side))
			graph.createBond(side, segment.getSide(ancestor.bond.left))
		else:	
			print 'Pulling up junction bond'
			# Pulling up
			liftedEdges = [X[0] for X in side.liftedPartners()]
			target = getMajority(liftedEdges)
			graph.validate()

			if graph.sideThread(side) is graph.sideThread(target):
				# Same thread => must connect side to target directly
				if target.bond is None:
					print 'Direct graft'
					graph.createBond(side, target)
					graph.validate()
				else:
					print 'Stubbing simulatenous bond confusion'
					graph.createBond(side, graph.newSegment().left)
					graph.validate()

			elif graph.eventGraph.testConstraint(graph.sideThread(target), graph.sideThread(side)):
				# if current side is not older than target create child
				# Note interpolating child could be a good idea, but we have coalescence correction to correct for that
				segment = graph.newSegment()
				graph.createBranch(target.segment, segment)
				graph.createBond(side, segment.getSide(target.left))
				graph.validate()

			else:
				# If target is definitely younger than current
				if target.parent() is None:
					# Aha you can give it a parent!
					segment = graph.newSegment()
					graph.createBranch(segment, target.segment)
					graph.createBond(side, target.parent())
					graph.validate()
				else:
					# Stubbing
					print 'Stubbing parent'
					graph.createBond(side, graph.newSegment().left)
					graph.validate()
			
	else:
		# Bridging
		children = filter(lambda X: X.bond is not None or len(X.liftedPartners()) > 0, side.children())
		if len(children) == 1:
			if side.parent() is None:
				graph.validate()
				print 'Root bridge'
				segment = graph.newSegment()
				graph.createBranch(side.segment, segment)
				graph.deleteBranch(side.segment, children[0].segment)
				graph.createBranch(segment, children[0].segment)

				segment2 = graph.newSegment()	
				graph.createBranch(side.bond.segment, segment2)
				graph.createBond(segment.getSide(side.left), segment2.getSide(side.bond.left))
				graph.validate()
			else:
				print 'Ignoring badness due to unattached junctions below'
			return
		children.pop(random.randrange(len(children)))
		if len(children) == 1 and children[0].bond is None and random.random() < 0.5:
			print 'Creating bridge bond on child'
			bridge = children[0]
		else:
			bridgeS = graph.newSegment()
			print 'Creating bridge bond on new segment', id(bridgeS)
			for child in children:
				graph.deleteBranch(side.segment, child.segment)
				graph.createBranch(bridgeS, child.segment) 
			graph.createBranch(side.segment, bridgeS)
			bridge = bridgeS.getSide(side.left)

		stepChildren = filter(lambda X: (X.bond is not None or len(X.liftedPartners()) > 0) and graph.eventGraph.testConstraint(graph.sideThread(bridge), graph.sideThread(X)), side.bond.children())
		if len(stepChildren) > 0:
			stepChild = random.choice(stepChildren)
			if graph.sideThread(stepChild) is graph.sideThread(bridge):
				if stepChild.bond is not None:
					print 'Stupid bridge stub'
					bridgePartnerS = graph.newSegment()
					graph.createBranch(side.bond.segment, bridgePartnerS)
				else:
					bridgePartnerS = stepChild.segment
			else:
				bridgePartnerS = graph.newSegment()
				graph.deleteBranch(side.bond.segment, stepChild.segment)
				graph.createBranch(bridgePartnerS, stepChild.segment) 
				graph.createBranch(side.segment, bridgePartnerS)
				print 'Interpolating bridge', id(bridgePartnerS)
		else:
			print 'Necessary bridge stub'
			bridgePartnerS = graph.newSegment()
			graph.createBranch(side.bond.segment, bridgePartnerS)
		graph.createBond(bridge, bridgePartnerS.getSide(side.bond.left))
		
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

detectors = [listCase1, listCase2, listCase3]
