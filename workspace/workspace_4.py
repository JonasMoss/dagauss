import networkx
from itertools import product

G = networkx.DiGraph()
G.add_nodes_from(["v", "w", "x", "y"])
G.add_edges_from([("v", "x"),
                  ("w", "y"),
                  ("x", "y")])
to_dag(G)

responses = ["x", "y"]
covariates = ["w"]
conditionants = ["v"]
conditionants = []
rsqsign(G, responses, covariates)

G = networkx.DiGraph()
G.add_nodes_from(["x", "y", "z", "w"])
G.add_edges_from([("x", "y"),
                  ("z", "x"),
                  ("w", "x")])
to_dag(G)

responses = ["y", "x"]
covariates = ["z", "w"]
conditionants = []

subs =  {"beta_x"  : 0,
         "beta_y"  : 0,
         "beta_z"  : 0,
         "beta_w"  : 0,
#         "beta_xy" : -1,
#         "beta_zx" : 1,
#         "beta_wx" : 1,
         "sigma_x" : 1,
         "sigma_y" : 1,
         "sigma_z" : 1,
         "sigma_w" : 1}

rsqsign(G, responses, covariates).subs(subs)

(sympy.simplify(A[:, 0].dot(B[:, 0]) + A[:, 1].dot(B[:, 1]))/(A.norm()*B.norm())).subs(subs)

beta(G, responses, covariates)
beta(G, covariates, responses)
covariance(G, responses, covariates)
covariance(G, covariates, responses)
rsquared(G, responses, covariates, conditionants).subs(subs)
rsquared(G, covariates, responses, conditionants).subs(subs)


A = sympy.Matrix([[sympy.cos("theta"), -sympy.sin("theta")], 
                  [sympy.sin("theta"), sympy.cos("theta")]])

sympy.trace(A*top*A.T)/sympy.trace(A*bottom*A.T)
sympy.trace(top)/sympy.trace(bottom)

(sympy.trace(top)/sympy.trace(bottom)).subs(subs)
(top.norm("fro")/bottom.norm("fro")).subs(subs)
(top.norm(2)/bottom.norm(2)).subs(subs)

rsq.subs(subs)
rsq2 = rsquared(G, ["x"], covariates, )
rsq2.subs(subs)


beta(G, responses, covariates).subs(subs)
covariance(G, responses, covariates).subs(subs)
beta(G,  ["x"], covariates).subs(subs)

quotient.subs(subs)

sympy.simplify(Rsq.subs(subs)/Rsq2.subs(subs))

responses = ["v", "w"]
covariates = ["x", "y"]
conditionants = []