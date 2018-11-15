import networkx as nx

# Add a DAG with confounder z.

G = nx.DiGraph()

focus_node = ["A"]
inblanket = ["I_1", "I_2", "I_3", "I_4", "I_5", "I_6"]
exblanket = ["O_1", "O_2", "O_3", "O_4", "O_5", "O_6"]

edges = [("O_1", "I_1"), ("O_2", "I_2"), ("I_1", "A"), ("I_2", "A"), 
         ("I_3", "O_3"), ("I_3", "I_4"), ("A", "I_4"), ("A", "I_5"), 
         ("I_5", "O_4"), ("I_5", "O_5"), ("I_6", "I_5"), ("O_6", "I_6")]

G.add_nodes_from(focus_node + exblanket + inblanket)
G.add_edges_from(edges)

# Calculate its non-DAG:
G = populate_attributes(G)
H = calculate_H(G)

# Find its variance and covariance.
mean(H)
covariance(H)

variables = focus_node 
conditionants = inblanket 

cond = conditionals(focus_node , inblanket, H)
