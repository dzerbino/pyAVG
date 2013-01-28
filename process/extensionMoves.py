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

####Make G-bounded AVG from G-bounded DNA history graph

def listCase1(graph):
	return [ExtensionMove(applyCase1, (segment, graph)) for segment in filter(lambda X: X.substitutionAmbiguity() > 0, graph.segments)]

def applyCase1(args):
	rootSegment, graph = args
	assert rootSegment.substitutionAmbiguity() > 0
	bottomSegment = random.choice(list(rootSegment.nonTrivialLiftedLabels()))
	#Get segments from chosen segment creating non-trivial label to root segment
	l = [ bottomSegment ]
	x = bottomSegment
	while x != rootSegment:
		x = x.parent
		l.append(x)
	#Now walk this list from the root and randomly label the first eligible segment
	while len(l) > 0:
		x = l.pop()
		if rootSegment.label == None or x == bottomSegment or len(rootSegment.liftedLabels()) - len(x.liftedLabels()) > 0:
			if x.label != None:
				assert x == bottomSegment
				print 'Adding necessary bridge label'
				graph.interpolateSegment(x).setLabel(str(rootSegment.label))
			elif len(x.liftedLabels()) == 1:
				print 'Adding necessary bridge label'
				x.setLabel(str(rootSegment.label))
			else:
				print 'Adding junction label'
				x.setLabel(str(random.choice(list(x.liftedLabels()) + [ rootSegment ]).label))
			break

###############################################
## Case 2
## Adding necessary bonds
###############################################

def createNecessaryBridge(rootSide, sideToAttach):
	pass

def applyCase2(args):
	rootSide, graph = args
	assert rootSide.rearrangementAmbiguity() > 0
	bottomSide = random.choice(list(rootSide.nonTrivialLiftedBonds()))
	assert bottomSide != rootSide
	#Get segments from chosen segment creating non-trivial label to root segment
	l = [ bottomSide ]
	x = bottomSide
	while x != rootSide:
		x = x.parent()
		l.append(x)
	#Now walk this list root and randomly label the first eligible segment
	while len(l) > 0:
		x = l.pop()
		if x.bond != None:
			assert x == bottomSide
			assert rootSide.bond != None
			print 'Adding necessary bridge bond'
			createNecessaryBridge(rootSide, graph.interpolateSegment(x.segment).getSide(x.left))
		elif len(x.liftedLabels()) > 1:
			print 'Adding junction label'
			l2 = []
			for y in x.liftedLabels():
				z = getConcomittantUnattachedSegment(y.bond, x)
				if z != None:
					l2.append(z)
			i = random.choice(xrange(len(l2) + 1))
			if i < len(l2):
				if l2[i].bond == None:
					x.createBond(x, l2[i])
				else:
					x.createBond(x, graph.interpolateSegment(l2[i].segment).getSide(l2[i].left))
			else:
				#We can either create a bridge like bond, or
				#we can attach to another segment that has non-trivial lifts
				#or to a segment that has no
				createNecessaryBridge(rootSide, x)
		else:
			assert x != rootSide
			print 'Adding necessary bridge bond'
			createNecessaryBridge(rootSide, x)

def listCase2(graph):
	return [ExtensionMove(applyCase2, (side, graph)) for side in filter(lambda X: X.rearrangementAmbiguity() > 0, graph.sides())]

def getFirstUnattachedJunctionSideInFace(module):
	pass

def getFirstAttachedSide(side):
	pass

def getPartnerJunctions(side):
	pass

def applyCase2(args):
	module, graph = args
	assert module.rearrangementAmbiguity() > 0
	x = getFirstUnattachedJunctionSideInFace(module)
	if x != None:
		print 'Attaching junction side'
		graph.createBond(x, random.choice(getPartnerJunctions(x)))
		#Remove any ping-pong bonds created..
	else:
		for side in module.sides:
			for child in side.children:
				x = getFirstAttachedSide(child)
				if x != None and x.bond.ancestor() != side.bond: 
					print 'Adding necessary bridge bond'
					graph.createBond(graph.interpolateSegment(x).getSide(x.left), )
					break
		else:
			raise RuntimeError("Did not perform extension for module with rearrangement ambiguity")

def applyCase2(args):
	side, graph = args
	print 'Resolving rearrangement ambiguity at', id(side.segment)
	if side.bond is None:
		# In an ideal world Fitch parsimony would be a nice touch...
		ancestor = side.ancestor()
		if ancestor is not None and ancestor.bond is not None:
			print 'Pulling down junction bond'
			segment = graph.newSegment()
			graph.createBranch(ancestor.bond.segment, segment)
			graph.createBond(side, segment.getSide(ancestor.bond.left))
		else:	
			print 'Pulling up junction bond'
			liftedEdges = side.liftedPartners()
			target = getMajority(X[0] for X in liftedEdges)

			if graph.sideThread(side) is graph.sideThread(target):
				# Same thread => must connect side to target directly
				if target.bond is None:
					print 'Direct graft'
					graph.createBond(side, target)
				else:
					print 'Stubbing simulatenous bond confusion'
					graph.createBond(side, graph.newSegment().left)

			elif target.bond is None and graph.eventGraph.testConstraint(graph.sideThread(target), graph.sideThread(side)) and graph.eventGraph.testConstraint(graph.sideThread(side), graph.sideThread(target)):
				print 'Direct graft, separate threads'
				graph.createBond(side, target)

			elif graph.eventGraph.testConstraint(graph.sideThread(target), graph.sideThread(side)):
				print 'Bridging along lifting path'
				# if current side is not older than target interpolate child
				# Start by looking for child whose bond lifts to target
				children = filter(lambda X: X[0] is target, liftedEdges)
				child = random.choice(children)[2]
				# Follow lifting path up towards target
				ptr = child.bond.parent()
				while not graph.eventGraph.testConstraint(graph.sideThread(ptr), graph.sideThread(side)):
					ptr = ptr.parent()
				# Bridge
				if graph.sideThread(ptr) is graph.sideThread(side):
					graph.createBond(side, ptr)
					print 'Found sibling on lifting path', id(ptr)
				else:
					segment = graph.newSegment()
					print 'Interpolated sibling on lifting path', id(segment)
					graph.createBranch(ptr.segment, segment)
					graph.createBond(side, segment.getSide(target.left))

			else:
				# If target is definitely younger than current
				if target.parent() is None:
					# Aha you can give it a parent!
					print 'Adding parentage'
					segment = graph.newSegment()
					graph.createBranch(segment, target.segment)
					graph.createBond(side, target.parent())
				else:
					# Stubbing
					print 'Stubbing parent'
					graph.createBond(side, graph.newSegment().left)
			
	else:
		# Bridging
		children = filter(lambda X: X.bond is not None or len(X.liftedPartners()) > 0, side.children())
		if len(children) == 1:
			if side.parent() is None:
				print 'Root bridge'
				segment = graph.newSegment()
				graph.createBranch(side.segment, segment)
				graph.deleteBranch(side.segment, children[0].segment)
				graph.createBranch(segment, children[0].segment)

				segment2 = graph.newSegment()	
				graph.createBranch(side.bond.segment, segment2)
				graph.createBond(segment.getSide(side.left), segment2.getSide(side.bond.left))
			else:
				print 'Ignoring badness due to unattached junctions below'
			return
		children.pop(random.randrange(len(children)))
		if len(children) == 1 and children[0].bond is None and random.random() < 0.9:
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
