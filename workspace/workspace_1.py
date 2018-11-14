import networkx as nx
import sympy

G = nx.DiGraph()
G.add_node("x", sigma = sympy.Symbol("sigma_x"), beta = sympy.Symbol("beta_x"))
G.add_node("y", sigma = sympy.Symbol("sigma_y"), beta = sympy.Symbol("beta_y"))
G.add_node("z", sigma = sympy.Symbol("sigma_z"), beta = sympy.Symbol("beta_z"))
G.add_node("u", sigma = sympy.Symbol("sigma_u"), beta = sympy.Symbol("beta_u"))
G.add_edge("x", "y", beta = sympy.Symbol("beta_xy"))
G.add_edge("x", "z", beta = sympy.Symbol("beta_xz"))
G.add_edge("u", "x", beta = sympy.Symbol("beta_ux"))
G.add_edge("u", "y", beta = sympy.Symbol("beta_uy"))
G.add_edge("y", "z", beta = sympy.Symbol("beta_yz"))


nx.is_directed_acyclic_graph(G)


G = nx.DiGraph()
G.add_node("x", sigma = sympy.Symbol("sigma_x"), beta = sympy.Symbol("beta_x"))
G.add_node("y", sigma = sympy.Symbol("sigma_y"), beta = sympy.Symbol("beta_y"))
G.add_node("z", sigma = sympy.Symbol("sigma_z"), beta = sympy.Symbol("beta_z"))
G.add_node("u", sigma = sympy.Symbol("sigma_u"), beta = sympy.Symbol("beta_u"))
G.add_edge("x", "y", beta = sympy.Symbol("beta_xz"))
G.add_edge("u", "x", beta = sympy.Symbol("beta_ux"))
G.add_edge("u", "z", beta = sympy.Symbol("beta_uy"))
G.add_edge("z", "y", beta = sympy.Symbol("beta_yz"))


G = nx.DiGraph()
G.add_node("x", sigma = sympy.Symbol("sigma_x"), beta = sympy.Symbol("beta_x"))
G.add_node("y", sigma = sympy.Symbol("sigma_y"), beta = sympy.Symbol("beta_y"))
G.add_node("z", sigma = sympy.Symbol("sigma_z"), beta = sympy.Symbol("beta_z"))
G.add_edge("x", "y", beta = sympy.Symbol("beta_xy"))
G.add_edge("z", "x", beta = sympy.Symbol("beta_zx"))
G.add_edge("z", "y", beta = sympy.Symbol("beta_zy"))

G = nx.DiGraph()
G.add_nodes_from(["x", "y", "z"])
G.add_edges_from([("x", "y"),
                  ("z", "x"),
                  ("z", "y")])

G = populate_attributes(G)
H = calculate_H(G)

mean(H)
covariance(H)

mean(H).subs(
        {"beta_z" : 1,
         "beta_y" : 1,
         "beta_x" : 1})

covariance(H).subs(
        {"beta_zy" : 1,
         "beta_zx" : 1,
         "sigma_z" : 1,
         "sigma_x" : 1,
         "sigma_y" : 1})


variables = ["y"]
conditionants = ["x", "z"]

cond = conditionals(variables, conditionants, H)
cond["mean"]
cond["cov"]

var = sympy.collect(sympy.expand(cond["mean"][0]), "x").coeff("x", 1)
sympy.cancel(var).subs(
        {"beta_z" : 1,
         "beta_y" : 1,
         "beta_x" : 1,
         "beta_zy" : 1,
         "beta_zx" : 1,
         "sigma_z" : 1,
         "sigma_x" : 1,
         "sigma_y" : 1})

