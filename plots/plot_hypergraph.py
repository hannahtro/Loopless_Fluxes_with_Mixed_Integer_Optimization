from halp.directed_hypergraph import DirectedHypergraph

# Initialize an empty hypergraph
H = DirectedHypergraph()

# Add nodes 's' and 't' individually with arbitrary attributes
H.add_node('s', source=True)
H.add_node('t', sink=True)
# Add several nodes simultaneously, having the same arbitrary attributes 
H.add_nodes(['x', 'y', 'z', 'u', 'a', 'b'], label='grey')

# Add hyperedge from {'s'} to {'x'} with a weight of 1
H.add_hyperedge(set(['s']), set(['x']), weight=1)
# Add hyperedge from {'s'} to {'x', 'y'} with some arbitrary attributes and weight of 2
H.add_hyperedge(set(['s']), set(['x', 'y']), {'color': 'red', 'active': True}, weight=2)
# Add several hyperedges simultaneously, having individual weights
hyperedges = [(['s'], ['z'], {'weight': 2}),
              (['s'], ['t'], {'weight': 100}),
              (['x'], ['s'], {'weight': 1}),
              (['x', 'y', 'z'], ['u', 't'], {'weight': 3}),
              (('t', 'b'), ('a'), {'weight': 1}),
              (set(['a']), set(['u', 't']), {'weight': 1})]
H.add_hyperedges(hyperedges)