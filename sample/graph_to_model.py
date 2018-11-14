# -*- coding: utf-8 -*-

from itertools import product
import sympy
import networkx 

sympy.init_printing(use_unicode = True)

# This should be made into a method of the DagNormal class.


def populate_attributes(G):
    """ Populates a directed graph G with attributes. """
    for node in G.nodes:
        G.nodes[node]["beta"] = sympy.Symbol("beta_" + node)
        G.nodes[node]["sigma"] = sympy.Symbol("sigma_" + node)
    for edge in G.edges:
        G.edges[edge]["beta"] = sympy.Symbol("beta_" + edge[0] + edge[1])
    return G

def calculate_H(G):
    """ Calculate the unconditional mean vector of a multivariate normal DAG.

    This function calculates the mean vector of a multivariate normal DAG. It
        only calculates the vector and stores the values in the dag. Use the
        function 'mean_vector' to get to the vector. 

    Args:
        G (DagNormal): A DagNormal object representing a multivariate normal.

    Returns:
        None: The function modifies G in place.

    Examples:
        Examples should be written in doctest format, and should illustrate how
        to use the function.

        >>> a = [1,2,3]
        >>> print [x + 3 for x in a]
        [4, 5, 6]

    """
    V = list(networkx.topological_sort(G))
    H = G.to_undirected()
    H.graph["sort"] = V

    # This loop takes care of the means of the unconditional model.
    for v in V:
        edges = list(G.in_edges(v))
        vertices = [x for (x, _) in edges]
        self_contribution = H.nodes[v]["beta"]
        parent_contribution = sum([H.nodes[vertex]["mu"]*H.edges[edge]["beta"] 
                                   for vertex, edge 
                                   in zip(vertices, edges)])
    
        H.nodes[v]["mu"] = self_contribution + parent_contribution

    # This loop takes care of sigmas of the unconditional model.
    for i, v in zip(range(0, len(V)), V):
        edges = list(G.in_edges(v))
        vertices = [x for (x, _) in edges]
        product_edges = product(vertices, vertices)
        parent_contribution = sum([H.edges[(x, y)]["psi"]*H.edges[(x, v)]["beta"]*H.edges[(y, v)]["beta"] 
                                  for (x, y) 
                                  in product_edges])
        self_contribution = H.nodes[v]["sigma"]**2
        H.add_edge(v, v, psi = self_contribution + parent_contribution)
        
        predecessors = V[:i]
        for w in predecessors:
            contribution = sum([H.edges[edge]["beta"]*H.edges[(edge[0], w)]["psi"] for edge in edges])
            H.add_edge(v, w, psi = contribution)
            
    return H


def mean(H):
    V = H.graph["sort"]
    return sympy.Matrix([H.nodes[v]["mu"] for v in V])

def covariance(H):
    V = H.graph["sort"]
    
    k = len(H)
    M = sympy.zeros(k, k)
    
    for (i, j) in product(range(0, k), range(0, k)):
        M[i, j] = H.edges[(V[i], V[j])]["psi"]
        
    return M
    
def conditionals(variables, conditionants, H):
    
    cov = covariance(H)
    mean_ = mean(H)
    
    V = H.graph["sort"]
    
    variables_indices = [V.index(var) for var in variables]
    conditionants_indices = [V.index(var) for var in conditionants]
    
    
    cov_AA = cov[variables_indices, variables_indices]
    cov_AB = cov[variables_indices, conditionants_indices]
    cov_BA = cov[conditionants_indices, variables_indices]
    cov_BB_inv = sympy.Inverse(cov[conditionants_indices, conditionants_indices])
    mean_A = sympy.Matrix([mean_[index] for index in variables_indices])
    mean_B = sympy.Matrix([mean_[index] for index in conditionants_indices])
    conditionants = sympy.Matrix(conditionants)
    
    # Fix ordering here!
    new_mean = mean_A + cov_AB*cov_BB_inv*(conditionants - mean_B)
    new_cov  = cov_AA - cov_AB*cov_BB_inv*cov_BA
    return {"mean": sympy.simplify(new_mean), 
            "cov": sympy.simplify(new_cov)}

def new_graph(H, G):
    """ Make a DAG with the same parameters the same parameters as an old one. 
    
    Raises an error of the conditional independence structure of H is 
        incompatible with that of G. 
        
    Args:
        H (networkx.DiGraph): An object inheriting from DiGraph fiving the 
        G (DagNormal): A DagNormal object representing a multivariate normal.
    """
    return None
    