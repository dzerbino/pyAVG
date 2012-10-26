#!/usr/bin/env python

class Module(object):
	def __init__(self, iter):
		self.sides = set(iter)
	
	def complexity_lb(self):
		pass	

	def complexity(self, lowerBound=True):
		if lowerBound:
			return self.complexity_lb()
		else:
			return self.complexity.ub()
