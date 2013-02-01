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
		l = rootSegment.liftedLabels()
		if rootSegment.label != None:
			l.add(rootSegment)
		x.setLabel(str(random.choice(list(l))))

###############################################
## Case 2
## Adding necessary bonds
###############################################

def sideAttachment(sidetoAttach, otherSideToAttach, graph):
	graph.createBond(sideToAttach, otherSideToAttach)

def branchAttachment(sideToAttach, branch, graph):
	graph.createBond(sideToAttach, graph.interpolateSegment(branch).getSide(branch.left))

def childAttachment(sideToAttach, parent, graph):
	x = graph.newSegment(graph)
	graph.createBranch(parent.segment, x)
	graph.createBond(sideToAttach, x.getSide(parent.left))
	
def stubAttachment(sideToAttach, nullArgument, graph):
	assert nullArgument == None
	graph.createBond(x, graph.newSegment().getSide(True))

def getUnattachedJunctionsAndEligibleUnattachedSidesOnLineage(bottomSide):
	#On the path from rootSide to bottomSide make list of junction sides
	#and sides with no unattached junction ancestors in this path
	l = []
	x = bottomSide.parent()
	#First get the eligible unattached junctions
	while x != None and x.bond == None:
		if len(x.liftedBonds()) > 1:
			l.append(x)
		x = x.parent()
	if x != None:
		assert x.bond != None #If root is unattached then it can be attached
		if len(l) > 0:
			x = l[-1].parent()
		else:
			x = bottomSide
		while x.bond == None or x == bottomSide:
			assert x == bottomSide or len(x.liftedBonds()) == 1
			l.append(x)
			x = x.parent()
	else:
		assert len(l) > 0
	return l

def getUnattachedSidesOnSideLineage(bottomSide):
	l = []
	while bottomSide != None and bottomSide.bond == None:
		l.append(bottomSide)
		bottomSide = bottomSide.parent()
	return l

def getConcomittantPartners(bottomSide, sideToAttach, graph, eligibleSidesFn):
	l = []
	for x in eligibleSidesFn(bottomSide):
		if graph.threadCmp(sideToAttach.segment.thread(), x.segment.thread()) == 0:
			l.append((sideAttachment, x))
		if graph.threadCmp(sideToAttach.segment.thread(), x.segment.thread()) <= 0 and \
		(x.parent() == None or graph.threadCmp(sideToAttach.segment.thread(), x.segment.parent.thread()) >= 0):
			l.append((branchAttachment, x))
	return l

def getPossibleBridgeBonds(rootSide, sideToAttach, graph, eligibleSidesFn):
	assert rootSide.bond != None
	return reduce(lambda x, y : x + y, [ [ (childAttachment, rootSide.bond)] ] + [ getConcomittantPartners(z, sideToAttach, graph, eligibleSidesFn) for z in rootSide.bond.nonTrivialLiftedBonds() ])
		
def getPossibleJunctionBonds(sideToAttach, graph):
	return reduce(lambda x, y : x + y, [ [] ] + [ getConcomittantPartners(z.bond, sideToAttach, graph, getUnattachedSidesOnSideLineage) for z in sideToAttach.liftedBonds() ])

def chooseRandomUnattachedSegmentOnSideLineage(bottomSide, graph):
	x = random.choice(getUnattachedJunctionsAndEligibleUnattachedSidesOnLineage(bottomSide))
	#If we've chosen the bottom side then it is attached already, so we need to interpolate, else we do it for the hell of it..
	if x == bottomSide or (len(x.liftedBonds()) == 1 and random.random() > 0.5):
		x = graph.interpolateSegment(x.segment).getSide(x.left)
	return x

def applyCase2(args):
	rootSide, graph = args
	assert rootSide.rearrangementAmbiguity() > 0
	
	bottomSide = random.choice(list(rootSide.nonTrivialLiftedBonds()))
	assert bottomSide != rootSide
	
	sideToAttach = chooseRandomUnattachedSegmentOnSideLineage(bottomSide, graph)
	assert sideToAttach.bond == None
	
	#Now proceed to attach chosen side
	if len(sideToAttach.liftedBonds()) > 1:
		l = getPossibleJunctionBonds(sideToAttach, graph) 
		if rootSide.bond != None:
			l += getPossibleBridgeBonds(rootSide, sideToAttach, graph, getUnattachedSidesOnSideLineage)
	else:
		assert rootSide.bond != None
		l = getPossibleBridgeBonds(rootSide, sideToAttach, graph, getUnattachedJunctionsAndEligibleUnattachedSidesOnLineage)
	if len(l) == 0:
		assert rootSide.bond == None
		l.append((stubAttachment, None))
	i = random.choice(l)
	i[0](sideToAttach, i[1], graph)

def listCase2(graph):
	return [ExtensionMove(applyCase2, (side, graph)) for side in filter(lambda X: X.rearrangementAmbiguity() > 0, graph.sides())]
		
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
