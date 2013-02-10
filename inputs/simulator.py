#!/usr/bin/env python


import sys
import random
import copy

from pyAVG.DNAHistoryGraph.evoHist import EvolutionaryHistory
from pyAVG.DNAHistoryGraph.thread import CircularSequenceThread
from pyAVG.DNAHistoryGraph.thread import CircularThread
from pyAVG.DNAHistoryGraph.label import Label

"""Produces random evolutionary histories"""

BRANCHPROB = 0.1
""" Probability that a branch point occurs at a given node """

class HistoryBranch(object):
	"""Branch in a history"""
        #########################################
        ## Basics
        #########################################
	def __init__(self):
		self.children = []

	def __str__(self):
		return "\n".join([self._label()] + map(str, self.children))

        #########################################
        ## Stats
        #########################################
	def _enumerate(self):
		return [self] + sum((X._enumerate() for X in self.children), [])

	def _cost(self):
		return self._operationCost() + sum(X._cost() for X in self.children)

	def _subs(self):
		return self._substitutionCost() + sum(X._subs() for X in self.children)

        #########################################
        ## GraphViz representation
        #########################################
	def _dot(self):
		""" GraphViz output """
		return "\n".join(["digraph G {" , "node [shape=rectangle]", self._dot2(), "}"])

	def _dot2(self):
		""" GraphViz output """
		return "\n".join([self._dotString()] + [X._dot2() for X in self.children])

	def _dotString(self):
		""" GraphViz output """
		return "\n".join([self._dotLabel()] + [self._dotEdge(X) for X in self.children])

	def _dotLabel(self):
		label = str(self.genome)
		return '%i [label="%s"]' % (id(self), label)

	def _dotEdge(self, child):
		return '%i -> %i [label="%s"]' % (id(self), id(child), child._dotBlurb())

#########################################
## Top Branch
#########################################

class InitialBranch(HistoryBranch):
	"""Root branch of an evolutionary history"""

	def __init__(self, length):
		super(InitialBranch, self).__init__()
		self.genome = ["".join(random.choice([ "A", "C", "T", "G" ]) for X in range(1,length+ 1))]

	def _label(self):
		return "INIT\n" + str(self.genome)

	def _operationCost(self):
		return 0

	def _substitutionCost(self):
		return 0

	def threads(self):
		T = [CircularSequenceThread(self.genome[0])]
		return sum([X.threads(T) for X in self.children], []) + [T]

#########################################
## Transformation branches
#########################################

class Operation(HistoryBranch):
	"""Non-root branch in an evolutionary history"""

	def __init__(self, parent):
		super(Operation, self).__init__()
		self.parent = parent
		parent.children.append(self)
		self.genome = self._product(parent.genome)

	def threads(self, parentThreads):
		T = self.modifyThreads(map(lambda X: X.childThread(), parentThreads))
		return sum([X.threads(T) for X in self.children], []) + [T]

class Identity(Operation):
	"""Nothing happens"""
	def __init__(self, parent):
		super(Identity, self).__init__(parent)

	def _product(self, genome):
		return genome

	def _label(self):
		return "ID"

	def _dotBlurb(self):
		return "ID"

	def _operationCost(self):
		return 0

	def _substitutionCost(self):
		return 0

	def modifyThreads(self, threads):
		return threads

def complement(base):
	if base == 'A':
		return 'T'
	elif base == 'T':
		return 'A'
	elif base == 'C':
		return 'G'
	elif base == 'G':
		return 'C'
	else:
		assert False, 'Unknown base for complementation: %s' % base

def revcomp(sequence):
	return "".join(map(complement, reversed(sequence)))

def revcompThread(sequence):
	return [X.opposite() for X in reversed(sequence)]

class DCJ(Operation):
	"""DCJ branch"""
	def __init__(self, parent, chrA, posA, chrB, posB, orientation):
		if chrA < chrB:
			self.chrA = chrA
			self.posA = posA
			self.chrB = chrB
			self.posB = posB
		elif chrA > chrB:
			self.chrA = chrB
			self.posA = posB
			self.chrB = chrA
			self.posB = posA
		else:
			self.chrA = chrA
			self.posA = min(posA, posB)
			self.chrB = chrB
			self.posB = max(posA, posB)
		self.orientation = orientation
		super(DCJ, self).__init__(parent)

	def _operationCost(self):
		return 1

	def _substitutionCost(self):
		return 0

	def _product(self, genome):
		newGenome = copy.copy(genome)
		if self.chrA != self.chrB:
			chrA = newGenome.pop(self.chrA)
			chrB = newGenome.pop(self.chrB - 1)
			if self.orientation:
				newGenome.append(chrA[:self.posA] + chrB + chrA[self.posA:])
			else:
				newGenome.append(chrA[:self.posA] + revcomp(chrB) + chrA[self.posA:])
		else:
			chrA = newGenome.pop(self.chrA)
			if self.orientation:
				newGenome.append(chrA[:self.posA] + chrA[self.posB:])
				assert len(newGenome[-1]) > 0
				newGenome.append(chrA[self.posA:self.posB])
				assert len(newGenome[-1]) > 0, self._label()
			else:
				newGenome.append(chrA[:self.posA] + revcomp(chrA[self.posA:self.posB]) + chrA[self.posB:])
				assert len(newGenome[-1]) > 0
		return newGenome

	def _label(self):
		return "DCJ\t%i\t%i\t%i\t%i\t%s\n%s" % (self.chrA, self.posA, self.chrB, self.posB, str(self.orientation), str(self.genome))

	def _dotBlurb(self):
		return "DCJ %i,%i,%i,%i" % (self.chrA, self.posA, self.chrB, self.posB)

	def modifyThreads(self, threads):
		if self.chrA != self.chrB:
			threadA = threads.pop(self.chrA)
			threadB = threads.pop(self.chrB - 1)
			if self.orientation:
				threads.append(CircularThread(threadA[:self.posA] + threadB[:] + threadA[self.posA:]))
			else:
				threads.append(CircularThread(threadA[:self.posA] + revcompThread(threadB) + threadA[self.posA:]))
		else:
			threadA = threads.pop(self.chrA)
			if self.orientation:
				threads.append(CircularThread(threadA[:self.posA] + threadA[self.posB:]))
				threads.append(CircularThread(threadA[self.posA:self.posB]))
			else:
				threads.append(CircularThread(threadA[:self.posA] + revcompThread(threadA[self.posA:self.posB]) + threadA[self.posB:]))
		return threads
		
class Mutation(Operation):
	"""Substituion branch"""
	def __init__(self, parent, chr, pos):
		self.chr = chr
		self.pos = pos 
		super(Mutation, self).__init__(parent)

	def _product(self, genome):
		return genome[:self.chr] + [genome[self.chr][:self.pos] + complement(genome[self.chr][self.pos]) + genome[self.chr][self.pos + 1:]] + genome[self.chr + 1:]

	def _mutationStr(self):
		return self.parent.genome[self.chr][self.pos] + ">" + self.genome[self.chr][self.pos]

	def _label(self):
		return "MUT\t%i\t%i\t%s\n%s" % (self.chr, self.pos, self._mutationStr(), str(self.genome))

	def _dotBlurb(self):
		return "MUT %i,%i,%s" % (self.chr,self.pos, self._mutationStr())

	def _operationCost(self):
		return 0

	def _substitutionCost(self):
		return 1

	def modifyThreads(self, threads):
		threads[self.chr][self.pos].segment.label = Label(threads[self.chr][self.pos].segment, threads[self.chr][self.pos].segment.label.complement())
		return threads

#########################################
## Evolutionary History
#########################################
class History(object):
	"""Evolutionary history"""

	def __init__(self, root):
		self.root = root

	def cost(self):
		return self.root._cost()

	def subs(self):
		return self.root._subs()

	def enumerate(self):
		return self.root._enumerate()

	def __str__(self):
		return str(self.root)

	def dot(self):
		""" GraphViz output """
		return self.root._dot()

	def avg(self):
		return EvolutionaryHistory(self.root.threads())

#########################################
## Random Evolutionary History
#########################################
def _addChildBranch(branch, choice, noDupes=False):
	if len(branch.genome) == 0:
		return
	elif choice < 0.5:
		chr = random.randrange(len(branch.genome))
		pos = random.randrange(len(branch.genome[chr]))
		Mutation(branch, chr, pos)
	elif choice < 0.9:
		chrA = random.randrange(len(branch.genome))
		posA = random.randrange(len(branch.genome[chrA]))
		chrB = random.randrange(len(branch.genome))
		posB = random.randrange(len(branch.genome[chrB]))
		if chrA == chrB and posA == posB:
			posB = (posA + 1) % len(branch.genome[chrA])
		orientation = (random.random() > 0.5)
		DCJ(branch, chrA, posA, chrB, posB, orientation)
	else:
		# Duplication
		Identity(branch)
		Identity(branch)

def _extendHistory(branch, counter, randomNums):
	_addChildBranch(branch, randomNums[counter-1])
	if counter > 1:
		map(lambda X: _extendHistory(X, counter - 1, randomNums), branch.children)

class RandomHistory(History):
	def __init__(self, length, maxDepth):
		randomNums = [random.random() for X in range(maxDepth)]
		root = InitialBranch(length)
		_extendHistory(root, maxDepth, randomNums)
		super(RandomHistory, self).__init__(root)

#########################################
## Unit test
#########################################
def test_main():
	history = RandomHistory(5,5)
	print history
	print history.dot()
	avg = history.avg()
	print avg.dot()
	assert avg.validate()
	assert avg.isAVG(), avg.dot()
	assert history.subs() == avg.lowerBoundSubstitutionCost(), "%s\t%s" % (history.subs(), avg.substitutionCost())
	assert history.cost() == avg.lowerBoundRearrangementCost(), "\n".join([str(history), avg.dot(), str(history.cost()), str(avg.rearrangementCost())])
	assert avg.lowerBoundSubstitutionCost() == avg.upperBoundSubstitutionCost()
	assert avg.lowerBoundRearrangementCost() == avg.upperBoundRearrangementCost()

if __name__ == '__main__':
	test_main()
