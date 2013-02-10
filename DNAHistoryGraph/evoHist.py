from graph import DNAHistoryGraph

def _layerRanks(layer):
	return " ".join(["{rank=same;"] + [str(id(T.segment)) for Th in layer for T in Th] + ["}"])

class EvolutionaryHistory(DNAHistoryGraph):
	def __init__(self, layers):
		self.layers = layers
		super(EvolutionaryHistory, self).__init__(traversal.segment for layer in layers for thread in layer for traversal in thread)

	def dot(self):
		ranks = map(_layerRanks, self.layers)
		return "\n".join(["digraph G {"] + ranks + [X.dot() for X in self.segments] + ["}"])
