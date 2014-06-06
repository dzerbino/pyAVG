pyAVG
=====

AVG sampler prototype

Implementation of the DNA history data structure. Currently, implementation of the definition:
- Structure
- Ambiguity evaluation
- Complexity evaluation (upper and lower)
- Random generation of AVGs and DNA history graphs (->unit testing)
- Graphviz output
- Sampling of AVG extensions thru Ben`s algos

To run:
- Set your PYTHONPATH to find the current directory.
- python process/reAVG.py
  - Creates a random History
  - Creates an AVG (pseudo) realization of the history
  - Removes random elements from the AVG to produce a DNA history graph
  - Tries to construct an AVG extension of the DNA history graph using the heuristic moves

