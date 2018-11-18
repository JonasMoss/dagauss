# -*- coding: utf-8 -*-

from itertools import product
import sympy
import networkx 
import copy

sympy.init_printing(use_unicode = True)

# This should be made into a method.

def get_order(G, values = [], keep = False):
    
    H = G.graph["dependency_graph"]
    
    order = copy.deepcopy(H.graph["sort"])
    
    if not values: return order
    
    if keep == False:
        for value in values: 
            order.remove(value)
            
    else:
        for value in set(order).difference(values): 
            order.remove(value)
            
    return order


def variable_indices(G, values, restrictions = [], sort = False):
    if(not restrictions):
        restrictions = get_order(G)
    order = get_order(G, restrictions, keep = True)
    indices = [order.index(value) for value in values]
    if sort: indices.sort()
    return indices

# This is the ___init___ method.
    
def to_dag(G):
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

    """ Populates a directed graph G with attributes. """
    for node in G.nodes:
        G.nodes[node]["beta"] = sympy.Symbol("beta_" + node, real = True )
        G.nodes[node]["sigma"] = sympy.Symbol("sigma_" + node, positive = True)
    for edge in G.edges:
        G.edges[edge]["beta"] = sympy.Symbol("beta_" + edge[0] + edge[1], real = True)

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
    G.graph["dependency_graph"] = H
    
# This is the parameters() method.
def parameters(G, responses = [], covariates = []):
    
    if(not responses): 
        responses = get_order(G, covariates)
        
    H = G.graph["dependency_graph"]
    order = get_order(G)
    mean_ = sympy.Matrix([H.nodes[value]["mu"] for value in order])
    cov = sympy.zeros(len(order), len(order))
    for (i, j) in product(range(len(order)), range(len(order))):
        cov[i, j] = H.edges[(order[i], order[j])]["psi"]

    if(not covariates):
        return {"mean": mean_[variable_indices(G, responses, sort = True), 0], 
                "cov": cov[variable_indices(G, responses, sort = True), 
                           variable_indices(G, responses, sort = True)]}
    
    V = get_order(G)
    
    responses_indices = variable_indices(G, responses, sort = True)
    covariates_indices = variable_indices(G, covariates, sort = True)
    
    cov_AA = cov[responses_indices, responses_indices]
    cov_AB = cov[responses_indices, covariates_indices]
    cov_BA = cov[covariates_indices, responses_indices]
    cov_BB_inv = sympy.Inverse(cov[covariates_indices, covariates_indices])
    mean_A = sympy.Matrix([mean_[index] for index in responses_indices])
    mean_B = sympy.Matrix([mean_[index] for index in covariates_indices])

    new_mean = mean_A + cov_AB*cov_BB_inv*(sympy.Matrix(covariates) - mean_B)
    new_cov  = cov_AA - cov_AB*cov_BB_inv*cov_BA
    return {"mean": sympy.simplify(new_mean), 
            "cov": sympy.simplify(new_cov)}

# This is the mean() method.
def mean(G, responses = [], covariates = []): 
    return parameters(G, responses = responses, covariates = covariates)["mean"]
    
# This is the covariance() method.
def covariance(G, responses = [], covariates = []): 
    return parameters(G, responses = responses, covariates = covariates)["cov"]

# The variance method picks the only item from the covariance matrix.
def variance(G, responses = [], covariates = []):
    if len(responses) == 1:
        return parameters(G, responses = responses, 
                             covariates = covariates)["cov"][0]
    else:
        return parameters(G, responses = responses, 
                             covariates = covariates)["cov"]
        
def beta(G, responses = [], covariates = [], conditionants = []):
    variables = covariates + conditionants    
    
    conditional_means = mean(G, responses, variables)
    
    def collect(index):
        expr = sympy.expand(conditional_means[index])
        return sympy.collect(expr = expr, 
                             syms = variables)
        
    collections = [collect(index) for index in range(len(conditional_means))]
    
    betas = sympy.Matrix([collection.coeff(variables) for 
                          collection, variables in 
                          product(collections, variables)])
    
    betas.reshape(len(collections), len(variables)).T
    
    indices = variable_indices(G, values = covariates, 
                                  restrictions = variables, 
                                  sort = True)
    
    return betas[indices, :]

def rsquared(G, responses, covariates, conditionants = [], norm = "trace"):
    """ Calculates the theoretical R squared. 

    This function calculates R squared, also known as the coefficient of 
       determination. 

    Args:
        G (DaGauss): A DaGauss object representing a multivariate normal.
        responses: The 

    Returns:
        None: The function modifies G in place.

    Examples:
        Examples should be written in doctest format, and should illustrate how
        to use the function.

        >>> a = [1,2,3]
        >>> print [x + 3 for x in a]
        [4, 5, 6]
    """
    
    cov = variance(G, responses = covariates, 
                      covariates = conditionants)

    betas = beta(G, responses = responses, 
                    covariates = covariates,
                    conditionants = conditionants)

    top = betas.T*cov*betas
    
    bottom = covariance(G, responses = responses, 
                           covariates = conditionants)
    
    if(norm == "trace"):
        return sympy.trace(top)/sympy.trace(bottom)
    else:
        return top.norm(norm)/bottom.norm(norm)
    


def correlation(G, variables = [], conditionants = []):
    """ Calculates the conditional correlation between variable1 and variable2
        given the conditionants. 

    """
    cov = covariance(G, variables, conditionants)
    k = cov.shape[0]
    sds = sympy.Matrix([1/sympy.sqrt(cov[i, i]) for i 
                        in range(0, k)]*k).reshape(k, k)
    
    cor = cov.multiply_elementwise(sds).multiply_elementwise(sds.T)
    return cor.applyfunc(sympy.simplify)

