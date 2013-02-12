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

def chooseRandomUnlabeledSegment(rootSegment, bottomSegment):
	#Get segments from chosen segment creating non-trivial label to root segment
	l = []
	x = bottomSegment
	while x != rootSegment:
		if rootSegment.label != None or (x != bottomSegment and len(x.liftedLabels()) > 1):
			l.append(x)
		x = x.parent
	if rootSegment.label == None:
		l.append(rootSegment)
	return random.choice(l)

def applyCase1(args):
	rootSegment, graph = args
	assert rootSegment.substitutionAmbiguity() > 0
	bottomSegment = random.choice(list(rootSegment.nonTrivialLiftedLabels()))
	x = chooseRandomUnlabeledSegment(rootSegment, bottomSegment)
	if x == bottomSegment:
		assert rootSegment.label != None
		assert x.label != None
		print 'Adding necessary bridge label'
		graph.interpolateSegment(x).setLabel(str(rootSegment.label))
	elif len(x.liftedLabels()) == 1:
		print 'Adding necessary bridge label'
		assert rootSegment.label != None
		x.setLabel(str(rootSegment.label))
	else:
		print 'Adding junction label'
		assert len(x.liftedLabels()) > 1
		l = rootSegment.liftedLabels()
		if rootSegment.label != None:
			l.add(rootSegment)
		labelToAdd = random.choice(list(l)).label
		assert labelToAdd != None
		x.setLabel(str(labelToAdd))

###############################################
## Case 2
## Adding necessary bonds
###############################################

def sideAttachment(sideToAttach, otherSideToAttach, graph):
	print "Doing side attachment"
	#assert graph.threadCmp(graph.sideThread(sideToAttach), graph.sideThread(otherSideToAttach)) == 0
	assert graph.areSiblings(graph.sideThread(sideToAttach), graph.sideThread(otherSideToAttach))
	graph.createBond(sideToAttach, otherSideToAttach)

def branchAttachment(sideToAttach, branch, graph):
	print "Doing branch attachment"
	graph.createBond(sideToAttach, graph.interpolateSegment(branch.segment).getSide(branch.left))

def childAttachment(sideToAttach, parent, graph):
	print "Doing child attachment"
	x = graph.newSegment()
	assert x.label == None
	graph.createBranch(parent.segment, x)
	graph.createBond(sideToAttach, x.getSide(parent.left))
	
def stubAttachment(sideToAttach, nullArgument, graph):
	print "Doing stub attachment"
	assert nullArgument == None
	graph.createBond(sideToAttach, graph.newSegment().getSide(True))

def getUnattachedJunctionSidesOnLineage(bottomSide):
	#On the path from rootSide to bottomSide make list of junction sides
	#and sides with no unattached junction ancestors in this path
	l = []
	x = bottomSide.parent()
	#First get the eligible unattached junctions
	while x != None and x.bond == None:
		if x.isJunction():
			l.append(x)
		x = x.parent()
	return l

def getUnattachedPotentialBridgeAndJunctionSidesOnLineage(bottomSide):
	l = getUnattachedJunctionSidesOnLineage(bottomSide)
	if len(l) > 0:
		x = l[-1].parent()
	else:
		x = bottomSide
	while x != None and (x.bond == None or x == bottomSide):
		l.append(x)
		x = x.parent()
	return l

def getAllUnattachedSidesOnLineage(bottomSide):
	l = []
	assert bottomSide.bond != None
	x = bottomSide.parent()
	while x != None and x.bond == None:
		l.append(x)
		x = x.parent()
	return l

def branchIsConcomittantWithSide(sideToAttach, branch, graph):
	return (branch.parent() == None or graph.eventGraph.testConstraint(graph.sideThread(branch.parent()), graph.sideThread(sideToAttach), True)) and \
		graph.eventGraph.testConstraint(graph.sideThread(sideToAttach), graph.sideThread(branch), True)

def getConcomittantPartners(bottomSide, sideToAttach, graph, eligibleSidesFn):
	l = []
	for x in eligibleSidesFn(bottomSide):
		if x.bond == None and graph.areSiblings(graph.sideThread(sideToAttach), graph.sideThread(x)): #graph.threadCmp(graph.sideThread(sideToAttach), graph.sideThread(x)) == 0:
			l.append((sideAttachment, x))
		if branchIsConcomittantWithSide(sideToAttach, x, graph):
			l.append((branchAttachment, x))
		
	return l

def getPossibleBondsFromAttachedAncestor(rootSide, sideToAttach, graph, eligibleSidesFn):
	if rootSide.bond == None:
		return []
	return reduce(lambda x, y : x + y, [ [ (childAttachment, rootSide.bond)] ] + [ getConcomittantPartners(z, sideToAttach, graph, eligibleSidesFn) for z in rootSide.bond.nonTrivialLiftedBonds() ])
		
def getPossibleBondsFromAttachedDescendants(sideToAttach, graph, eligibleSidesFn):
	return reduce(lambda x, y : x + y, [ [] ] + [ getConcomittantPartners(z.bond, sideToAttach, graph, eligibleSidesFn) for z in sideToAttach.liftedBonds() ])

def chooseRandomUnattachedSegmentOnSideLineage(bottomSide, graph):
	x = random.choice(getUnattachedPotentialBridgeAndJunctionSidesOnLineage(bottomSide))
	#If we've chosen the bottom side then it is attached already, so we need to interpolate, else we do it for the hell of it..
	if x == bottomSide or (len(x.liftedBonds()) == 1 and random.random() > 0.5):
		print "interpolating segment"
		x = graph.interpolateSegment(x.segment).getSide(x.left)
	return x

def addPossibleRandomBonds(sideToAttach, graph, l):
	#This ensures we consider the possibility of every possible branch
	while random.random() > 0.5 or len(l) == 0:
		i = random.choice(xrange(len(graph.segments)+1))
		if i == len(graph.segments): #Ignore ping-pongs
			l.append((stubAttachment, None))
		else:
			continue
			x = random.choice(graph.segments)
			x = random.choice([ x.left, x.right])
			if x.parent() == None or len((set([x]) | x.liftedBonds()) & x.ancestor().nonTrivialLiftedBonds()) > 0: #Avoid breaking a window
				if x.bond != None or random.random() > 0.5:
					if branchIsConcomittantWithSide(sideToAttach, x, graph):
						l.append((branchAttachment, x))
				elif graph.areSiblings(graph.sideThread(sideToAttach), graph.sideThread(x)):
					l.append((sideAttachment, x))
	return l

def applyCase2(args):
	rootSide, graph = args
	assert rootSide.rearrangementAmbiguity() > 0
	
	bottomSide = random.choice(list(rootSide.nonTrivialLiftedBonds()))
	assert bottomSide != rootSide
	
	sideToAttach = chooseRandomUnattachedSegmentOnSideLineage(bottomSide, graph)
	assert sideToAttach.bond == None
	
	#Now proceed to attach chosen side
	if len(sideToAttach.liftedBonds()) > 1:
		l = getPossibleBondsFromAttachedDescendants(sideToAttach, graph, getAllUnattachedSidesOnLineage) + getPossibleBondsFromAttachedAncestor(rootSide, sideToAttach, graph, getAllUnattachedSidesOnLineage)
		addPossibleRandomBonds(sideToAttach, graph, l)
	else:
		l = getPossibleBondsFromAttachedDescendants(sideToAttach, graph, getUnattachedJunctionSidesOnLineage) + getPossibleBondsFromAttachedAncestor(rootSide, sideToAttach, graph, getUnattachedPotentialBridgeAndJunctionSidesOnLineage)
	
	i = random.choice(l)
	i[0](sideToAttach, i[1], graph)

def listCase2(graph):
	return [ExtensionMove(applyCase2, (side, graph)) for side in filter(lambda X: X.rearrangementAmbiguity() > 0, graph.sides())]
		
###############################################
## Case 3
## Adding necessary coalescences
###############################################

def listCase3(graph):
	return [ExtensionMove(applyCase3, (segment, graph)) for segment in filter(lambda X: X.ambiguity() > 0 and len(X.children) > 2, graph.segments)]

def applyCase3(args):
	segment, graph = args 
	children = list(segment.children)
	print "Creating pull down"
	children.pop(random.randrange(len(children)))	
	bridge = graph.newSegment()
	for child in children:
		graph.deleteBranch(segment, child)
		graph.createBranch(bridge, child) 
	graph.createBranch(segment, bridge)
	#Needs to check that it hasn't created any non-minimal bonds

detectors = [listCase1, listCase2, listCase3]
