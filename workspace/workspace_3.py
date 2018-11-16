import networkx
from itertools import product

G = networkx.DiGraph()
G.add_nodes_from(["v", "w", "x", "y", "z"])
G.add_edges_from([("w", "x"),
                  ("w", "y"),
                  ("z", "x"),
                  ("z", "w"),
                  ("v", "w"),
                  ("v", "x")])

G = populate_attributes(G)
H = calculate_H(G)

responses = ["x", "y"]
covariates = ["w", "z"]
conditionants = ["v"]
#rsquared(H, responses, covariates)

cond_Rsq = rsquared(H, responses, covariates, conditionants)
uncond_Rsq = rsquared(H, responses, covariates, [])

subs =  {"beta_v"  : 0,
         "beta_w"  : 0,
         "beta_x"  : 0,
         "beta_y"  : 0,
         "beta_z"  : 0,
         "beta_wx" : 1,
         "beta_wy" : 1,
         "beta_zx" : 1,
         "beta_zw" : 1,
         "beta_xy" : 1,
         "beta_zy" : 1,
         "beta_vw" : 1,
         "sigma_v" : 1,
         "sigma_w" : 1,
         "sigma_x" : 1,
         "sigma_y" : 1,
         "sigma_z" : 1}

cond_Rsq.subs(subs)/uncond_Rsq.subs(subs)


[sympy.N(cond_Rsq.subs(subs)), sympy.N(uncond_Rsq.subs(subs))]
