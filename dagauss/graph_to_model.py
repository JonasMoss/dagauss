# -*- coding: utf-8 -*-

from itertools import product
import sympy
import networkx 

sympy.init_printing(use_unicode = True)

# This should be made into a method of the DagNormal class.


def populate_attributes(G):
    """ Populates a directed graph G with attributes. """
    for node in G.nodes:
        G.nodes[node]["beta"] = sympy.Symbol("beta_" + node, real = True )
        G.nodes[node]["sigma"] = sympy.Symbol("sigma_" + node, positive = True)
    for edge in G.edges:
        G.edges[edge]["beta"] = sympy.Symbol("beta_" + edge[0] + edge[1], real = True)
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
            
    H.graph["sort"] = list(set(V))
    H.graph["sort"].sort()
    return H



    
def cov_and_mean(H, responses = [], covariates = []):
    
    def unconditional_covariance(H, responses = []):
        if(not responses): responses = H.graph["sort"]
        V = list(set(H.graph["sort"]).intersection(set(responses)))
    
        M = sympy.zeros(len(V), len(V))
    
        for (i, j) in product(range(len(V)), range(len(V))):
            M[i, j] = H.edges[(V[i], V[j])]["psi"]
        
        return M

    def unconditional_mean(H, responses = []):
        if(not responses): responses = H.graph["sort"]
        V = set(H.graph["sort"]).intersection(set(responses))
        return sympy.Matrix([H.nodes[v]["mu"] for v in V])
         
    if(not responses): 
        responses = list(set(H.graph["sort"]).difference(set(covariates)))
    
    if(not covariates):
        return {"mean": unconditional_mean(H, responses), 
                "cov": unconditional_covariance(H, responses)}
        
    cov = unconditional_covariance(H)
    mean_ = mean(H)
    V = H.graph["sort"]
    
    responses_indices = [V.index(var) for var in responses]
    covariates_indices = [V.index(var) for var in covariates]
    
    cov_AA = cov[responses_indices, responses_indices]
    cov_AB = cov[responses_indices, covariates_indices]
    cov_BA = cov[covariates_indices, responses_indices]
    cov_BB_inv = sympy.Inverse(cov[covariates_indices, covariates_indices])
    mean_A = sympy.Matrix([mean_[index] for index in responses_indices])
    mean_B = sympy.Matrix([mean_[index] for index in covariates_indices])
    covariates = sympy.Matrix(covariates)

    new_mean = mean_A + cov_AB*cov_BB_inv*(covariates - mean_B)
    new_cov  = cov_AA - cov_AB*cov_BB_inv*cov_BA
    return {"mean": sympy.simplify(new_mean), 
            "cov": sympy.simplify(new_cov)}


def mean(H, responses = [], covariates = []): 
    return cov_and_mean(H, responses = responses, covariates = covariates)["mean"]
    
def covariance(H, responses = [], covariates = []): 
    return cov_and_mean(H, responses = responses, covariates = covariates)["cov"]
