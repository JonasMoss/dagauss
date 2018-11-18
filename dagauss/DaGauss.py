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
def parameters(G, variables = [], conditionants = []):
    """ Calculate the conditional mean vector and covariance matrix

    Args:
        G: A DaGauss object representing a multivariate normal.
        variables: The variables in the regression. Can be more than one.
        conditionants: The variables to condition on.

    Returns:
        A dictionary containing the theoretical conditional mean vector and
        the theoretical conditional covariance matrix.
        
    """
    if(not variables): 
        variables = get_order(G, conditionants)
        
    H = G.graph["dependency_graph"]
    order = get_order(G)
    mean_ = sympy.Matrix([H.nodes[value]["mu"] for value in order])
    cov = sympy.zeros(len(order), len(order))
    for (i, j) in product(range(len(order)), range(len(order))):
        cov[i, j] = H.edges[(order[i], order[j])]["psi"]

    if(not conditionants):
        return {"mean": mean_[variable_indices(G, variables, sort = True), 0], 
                "cov": cov[variable_indices(G, variables, sort = True), 
                           variable_indices(G, variables, sort = True)]}
    
    V = get_order(G)
    
    variables_indices = variable_indices(G, variables, sort = True)
    conditionants_indices = variable_indices(G, conditionants, sort = True)
    
    cov_AA = cov[variables_indices, variables_indices]
    cov_AB = cov[variables_indices, conditionants_indices]
    cov_BA = cov[conditionants_indices, variables_indices]
    cov_BB_inv = sympy.Inverse(cov[conditionants_indices, conditionants_indices])
    mean_A = sympy.Matrix([mean_[index] for index in variables_indices])
    mean_B = sympy.Matrix([mean_[index] for index in conditionants_indices])

    new_mean = mean_A + cov_AB*cov_BB_inv*(sympy.Matrix(conditionants) - mean_B)
    new_cov  = cov_AA - cov_AB*cov_BB_inv*cov_BA
    return {"mean": sympy.simplify(new_mean), 
            "cov": sympy.simplify(new_cov)}

# This is the mean() method.
def mean(G, variables = [], conditionants = []):
    """ Calculate the conditional mean vector

    Args:
        G: A DaGauss object representing a multivariate normal.
        variables: The variables in the regression. Can be more than one.
        conditionants: The variables to condition on.

    Returns:
        The theoretical conditional mean vector
        
    """
    return parameters(G, variables = variables, conditionants = conditionants)["mean"]
    
# This is the covariance() method.
def covariance(G, variables = [], conditionants = []): 
    """ Calculate the conditional covariance matrix 

    Args:
        G: A DaGauss object representing a multivariate normal.
        variables: The variables in the regression. Can be more than one.
        conditionants: The variables to condition on.

    Returns:
        The theoretical conditional covariance matrix
        
    """
    return parameters(G, variables = variables, conditionants  = conditionants )["cov"]

# The variance method picks the only item from the covariance matrix.
def variance(G, variables = [], conditionants = []):
    """ Calculate the conditional covariance matrix 

    Args:
        G: A DaGauss object representing a multivariate normal.
        variables: The variables in the regression. Can be more than one.
        conditionants: The variables to condition on.

    Returns:
        The theoretical regression coefficient.
        
    """
    if len(variables) == 1:
        return parameters(G, variables = variables, 
                             conditionants = conditionants)["cov"][0]
    else:
        return parameters(G, variables = variables, 
                             conditionants = conditionants)["cov"]
        
def beta(G, responses = [], covariates = [], conditionants = []):
    """ Calculate the theoretical beta coefficient of a regression

    Args:
        G: A DaGauss object representing a multivariate normal.
        responses: The responses in the regression. Can be more than one.
        covariates: The covariates of the regression. 
        conditionants: The variables the regression is conditioned on.

    Returns:
        The theoretical regression coefficient.
        
    """
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
        G: A DaGauss object representing a multivariate normal.
        responses: The responses in the regression. Can be more than one.
        covariates: The covariates of the regression. 
        conditionants: The variables the regression is conditioned on.
        norm: Optional covariance matrix norm. Defaults to "trace", which is
        recommended.

    Returns:
        The caclulated R squared. A scalar sympy object.
        
    """
    
    cov = variance(G, variables = covariates, 
                      conditionants = conditionants)

    betas = beta(G, responses = responses, 
                    covariates = covariates,
                    conditionants = conditionants)

    top = betas.T*cov*betas
    
    bottom = covariance(G, variables = responses, 
                           conditionants = conditionants)
    
    if(norm == "trace"):
        return sympy.trace(top)/sympy.trace(bottom)
    else:
        return top.norm(norm)/bottom.norm(norm)
    


def correlation(G, variables = [], conditionants = []):
    """ Calculates the conditional correlation.
    
    
    Args:
        G: A DaGauss object representing a multivariate normal.
        variables: The variables you wish to find the correlation matrix for.
        conditionants: The variables the correlation matrix is conditioned on.

    Returns:
        A correlation matrix.

    """
    
    cov = covariance(G, variables = variables, 
                        conditionants = conditionants)
    k = cov.shape[0]
    sds = sympy.Matrix([1/sympy.sqrt(cov[i, i]) for i 
                        in range(0, k)]*k).reshape(k, k)
    
    cor = cov.multiply_elementwise(sds).multiply_elementwise(sds.T)
    return cor.applyfunc(sympy.simplify)

