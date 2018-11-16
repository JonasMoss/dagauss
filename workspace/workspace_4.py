import networkx
from itertools import product

G = networkx.DiGraph()
G.add_nodes_from(["v", "w", "x", "y"])
G.add_edges_from([("v", "x"),
                  ("w", "y"),
                  ("x", "y")])

G = populate_attributes(G)
H = calculate_H(G)

responses = ["x", "y"]
covariates = ["v", "w"]

Rsq = rsquared(H, responses, covariates, [])

subs =  {"beta_v"  : 0,
         "beta_w"  : 0,
         "beta_x"  : 0,
         "beta_y"  : 0,
         "beta_vx" : 1,
         "beta_wy" : 1,
         "beta_xy" : 1,
         "sigma_v" : 1,
         "sigma_w" : 1,
         "sigma_x" : 1,
         "sigma_y" : 1}

Rsq.subs(subs)
