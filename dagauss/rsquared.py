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
    

def rsqsign(G, responses, covariates, conditionants = []):
    cov = variance(G, responses = covariates, 
                       covariates = conditionants)

    indices = variable_indices(G, values = covariates, 
                              restrictions = covariates + conditionants, 
                              sort = True)
    
    betas = beta(G, responses = responses, 
                    covariates = covariates + conditionants)[indices, :]

    bottom = variance(G, responses = responses, 
                         covariates = conditionants + covariates)

    product_matrix = bottom*(cov*betas).T
    grand_sum = 0
    
    for i in product_matrix: grand_sum += i
    
    return grand_sum
