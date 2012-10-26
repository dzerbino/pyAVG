#!/usr/bin/env python

import sys
import random
import cnavg.avg.graph

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

	def _testPosition(self, position):
		if len(self.children):
			return any(child._testPosition(position) for child in self.children)
		else:
			return True

        #########################################
        ## GraphViz representation
        #########################################
	def _dot(self, weights):
		""" GraphViz output """
		return "\n".join(["digraph G {" , "node [shape=rectangle]", self._dot2(weights), "}"])

	def _dot2(self, weights):
		""" GraphViz output """
		return "\n".join([self._dotString(weights)] + [X._dot2(weights) for X in self.children])

	def _dotString(self, weights):
		""" GraphViz output """
		return "\n".join([self._dotLabel(weights)] + [self._dotEdge(X) for X in self.children])

	def _dotLabel(self, weights):
		label = str(self.genome)
		if len(self.children) == 0:
			label += " (%f)" % weights[self]
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
		self.genome = range(1,length+ 1)

	def _label(self):
		return "INIT\n" + str(self.genome)

	def _operationCost(self):
		return 0

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

	def _operationCost(self):
		if self._persists():
			return 1
		else:
			return 0

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

	def _persists(self):
		return True

class Inversion(Operation):
	"""Inversion branch"""
	def __init__(self, parent, start, length):
		self.start = start
		self.length = length
		super(Inversion, self).__init__(parent)

	def _product(self, genome):
		if self.start + self.length < len(genome):
			return genome[:self.start] + [-X for X in reversed(genome[self.start:self.start+self.length])] + genome[self.start + self.length:]
		else:
			return genome[:self.start + self.length - len(genome)] + [-X for X in reversed(genome[self.start + self.length - len(genome):self.start])] + genome[self.start:]

	def _label(self):
		return "INV\t%i\t%i\n%s" % (self.start, self.length, str(self.genome))

	def _dotBlurb(self):
		str = "INV %i,%i" % (self.start, self.length)
		if self._persists():
			return str
		else:
			return "*" + str

	def _testPosition(self, position):
		position = position % len(self.genome)
		if len(self.children):
			if position >= self.start and position < self.start + self.length:
				newPosition = 2*self.start + self.length - position - 1
			else:
				newPosition = position
			return any(child._testPosition(newPosition) for child in self.children)
		else:
			return True

	def _persists(self):
		return (self._testPosition(self.start) and self._testPosition(self.start + self.length)) or (self._testPosition(self.start - 1) and self._testPosition(self.start + self.length - 1))

class Duplication(Operation):
	"""Duplication branch"""
	def __init__(self, parent, start, length):
		self.start = start
		self.length = length
		super(Duplication, self).__init__(parent)

	def _product(self, genome):
		if self.start + self.length < len(genome):
			return genome[:self.start] + genome[self.start:self.start+self.length] + genome[self.start:]
		else:
			return genome[:self.start] + genome[self.start:] + genome[:self.start+self.length-len(genome)] + genome[self.start:]

	def _label(self):
		return "DUP\t%i\t%i\n%s" % (self.start, self.length, str(self.genome))

	def _dotBlurb(self):
		str = "DUP %i,%i" % (self.start, self.length)
		if self._persists():
			return str
		else:
			return "*" + str

	def _testPosition(self, position):
		position = position % len(self.genome)
		if len(self.children):
			if position >= self.start and position < self.start + self.length:
				return any(child._testPosition(position) for child in self.children) or any(child._testPosition(position + self.length) for child in self.children)
			if position >= self.start + self.length:
				return any(child._testPosition(position + self.length) for child in self.children)
			else:
				return any(child._testPosition(position) for child in self.children)
		else:
			return True

	def _persists(self):
		return any(self._testPosition(X) for X in range(self.start, self.start + self.length))

class Deletion(Operation):
	"""Deletion branch"""
	def __init__(self, parent, start, length):
		self.start = start
		self.length = length
		super(Deletion, self).__init__(parent)

	def _product(self, genome):
		if self.start + self.length < len(genome):
			return genome[:self.start] + genome[self.start+self.length:]
		else:
			return genome[(self.start + self.length) - len(genome):self.start]

	def _label(self):
		return "DEL\t%i\t%i\n%s" % (self.start, self.length, str(self.genome))

	def _dotBlurb(self):
		str = "DEL %i,%i" % (self.start, self.length)
		if self._persists():
			return str
		else:
			return "*" + str

	def _testPosition(self, position):
		if len(self.genome) == 0:
			return False
		position = position % len(self.genome)
		if position < self.start or position >= self.start + self.length:
			if len(self.children):
				if position >= self.start + self.length:
					return any(child._testPosition(position - self.length) for child in self.children)
				else:
					return any(child._testPosition(position) for child in self.children)
			else:
				return True
		else:
			return False

	def _persists(self):
		return self._testPosition(self.start - 1) or self._testPosition(self.start + self.length)

#########################################
## Evolutionary History
#########################################
class History(object):
	"""Evolutionary history"""

	def __init__(self, root):
		self.root = root

	def cost(self):
		return self.root._cost()

	def enumerate(self):
		return self.root._enumerate()

	def __str__(self):
		return str(self.root)

	def dot(self):
		""" GraphViz output """
		return self.root._dot(self.weights)

#########################################
## Weighted Evolutionary History
#########################################

def _removeInitialPath(avg):
	nodes = sorted(avg.nodes())
	previousNode = nodes[-1]
		
	for index in range(len(nodes)/2):
		nextNode = nodes[2 * index]
		avg.setLiftedEdge(previousNode, nextNode, -1)
		avg.changeSegment(nextNode, 0, -1)
		previousNode = avg[nextNode].twin
	return avg

class WeightedHistory(History):
	"""Evolutionary history with weighted branches"""
	def __init__(self, root, weights):
		super(WeightedHistory, self).__init__(root)
		self.weights = weights

	def _avg_Branch(self, avg, branch):
		if len(branch.children) == 0 and len(branch.genome) > 0:
			weight = self.weights[branch]
			nodes = sorted(avg.nodes())
			if branch.genome[-1] > 0:
				previousNode = nodes[2 * branch.genome[-1] - 1]
			else:
				previousNode = nodes[2 * -branch.genome[-1] - 2]
				
			for index in range(len(branch.genome)):
				if branch.genome[index] > 0:
					nextNode = nodes[2 * branch.genome[index] - 2]
				else:
					nextNode = nodes[2 * -branch.genome[index] - 1]
				
				avg.addLiftedEdge(previousNode, nextNode, weight)
				avg.changeSegment(nextNode, 0, weight)
				previousNode = avg[nextNode].twin
		return avg

	def avg(self):
		"""Produce resultant sequence graph with flow values"""
		avg = _removeInitialPath(cnavg.avg.graph.linearGraph(len(self.root.genome) * 2))
		return reduce(self._avg_Branch, self.enumerate(), avg)

#########################################
## Random Evolutionary History
#########################################
def _addChildBranch(branch):
	choice = random.random()
	start = random.randrange(len(branch.genome))

	if len(branch.genome) > 1 and choice < 0.7:
		length = random.randrange(1, len(branch.genome))
		Inversion(branch, start, length)
	elif len(branch.genome) > 1 and choice < 0.8:
		# Separate length prob for different operations (e.g. long distance duplications followed by deletions make things moot)
		length = int(random.expovariate(0.1))
		if length >= len(branch.genome):
			length = len(branch.genome) - 1
		if length == 0:
			length = 1
		Duplication(branch, start, length)
	elif len(branch.genome) > 1 and choice < 0.9:
		length = int(random.expovariate(0.1))
		if length >= len(branch.genome):
			length = len(branch.genome) - 1
		if length == 0:
			length = 1
		Deletion(branch, start, length)
	else:
		Identity(branch)

def _birthDeathModel():
	choice = random.random()
	if choice < BRANCHPROB:
		return 2
	else:
		return 1

def _extendHistory_Branch(branch):
	if len(branch.genome) == 0:
		return
	for i in range(_birthDeathModel()):
		_addChildBranch(branch)

def _extendHistory(branch, counter):
	newBranches = _extendHistory_Branch(branch)
	if counter > 1:
		map(lambda X: _extendHistory(X, counter - 1), branch.children)

class RandomHistory(WeightedHistory):
	def __init__(self, length, maxDepth):
		root = InitialBranch(length)
		_extendHistory(root, maxDepth)
		super(RandomHistory, self).__init__(root, None)

#########################################
## Random Weighted Evolutionary History
#########################################

class RandomWeightedHistory(RandomHistory):
	def _weightBranchesAtRandom(self):
		#return dict((X, random.randrange(1,1e6)) for X in self.enumerate() if len(X.children) == 0)
		return dict((X,1) for X in self.enumerate() if len(X.children) == 0)

	def _normalizeBranches(self, weights):
		total = float(sum(weights.values()))
		return dict((X, weights[X]/total) for X in weights)

	def _weightBranches(self):
		return self._normalizeBranches(self._weightBranchesAtRandom())

	def __init__(self, length, maxDepth):
		super(RandomWeightedHistory, self).__init__(length, maxDepth)
		# Ugly side effect needed to get around linear chain of inheritance
		self.weights = self._weightBranches()

#########################################
## Unit test
#########################################
def main():
	history = RandomWeightedHistory(5, 3)
	avg = history.avg()
	print history
	print history.cost()
	print avg
	print history.dot()

if __name__ == '__main__':
	main()
