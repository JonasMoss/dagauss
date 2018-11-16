def collect_betas(conditional_means, covariates):
    def collect(index):
        expr = sympy.expand(conditional_means[index])
        return sympy.collect(expr = expr, 
                             syms = covariates)
        
    collections = [collect(index) for index in range(len(conditional_means))]
    
    betas = sympy.Matrix([collection.coeff(covariates) for 
                          collection, covariates in 
                          product(collections, covariates)])
        
    return betas.reshape(len(collections), len(covariates)).T

def rsquared(H, responses, covariates, conditionants):
    """ Calculates the theoretical R squared. 

    This function calculates R squared, also known as the coefficient of 
       determination. 

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
    
    order = H.graph["sort"]
    
    responses = list(set(responses))
    conditionants = list(set(conditionants))
    covariates = list(set(covariates).difference(set(conditionants)))
    covariates_and_conditionants = list(set(covariates).union(set(conditionants)))
    
    cov_given_conditionants = covariance(H, responses = responses, covariates = conditionants)
    cov_covariates_given_conditionants = covariance(H, responses = covariates, covariates = conditionants)
    mean_given_covariates_and_conditionants = mean(H, responses, covariates_and_conditionants)
    betas = collect_betas(mean_given_covariates_and_conditionants, covariates)
    
    if(len(cov_covariates_given_conditionants) is 1):
        cov_expected_responses_given_convariates_and_conditionants = betas*betas.T*cov_covariates_given_conditionants[0]
    else:
        cov_expected_responses_given_convariates_and_conditionants  = betas.T*cov_covariates_given_conditionants*betas

    return sympy.trace(cov_expected_responses_given_convariates_and_conditionants)/sympy.trace(cov_given_conditionants )
    




def correlation(H, variables = [], conditionants = []):
    """ Calculates the conditional correlation between variable1 and variable2
        given the conditionants. 

    """
    
    if(not conditionants):
        cov = covariance(H, variables)
    else:
        cov = conditionals(variables, conditionants, H)["cov"]
        
    cov = cov*sympy.fraction(cov[0, 0])[1]
    k = cov.shape[0]
    sds = sympy.Matrix([1/sympy.sqrt(cov[i, i]) for i 
                        in range(0, k)]*k).reshape(k, k)
    
    cor = cov.multiply_elementwise(sds).multiply_elementwise(sds.T)
    return cor.applyfunc(sympy.simplify)
    
variables = ["y", "x"]
conditionants = ["z"]
condcor = conditional_correlation(H, variable1, variable2, conditionants)
condcor.subs( 
        {"beta_z"  : 0,
         "beta_y"  : 0,
         "beta_x"  : 0,
         "beta_yz" : 1,
         "sigma_z" : 1,
         "sigma_x" : 1,
         "sigma_y" : 1})


# With collider

collider = networkx.DiGraph()
collider.add_nodes_from(["x", "y", "z"])
collider.add_edges_from([("x", "y"),
                  ("x", "z"),
                  ("y", "z")])

collider = populate_attributes(collider)
collider_H = calculate_H(collider)

variables = ["y"]
conditionants = ["x", "z"]
collider_rsq_both = rsquared(collider_H, variables, conditionants).subs( 
        {"beta_z" : 0,
         "beta_y" : 0,
         "beta_x" : 0,
         "beta_yz" : 1,
         "sigma_z" : 1,
         "sigma_x" : 1,
         "sigma_y" : 1})
    
variables = ["y"]
conditionants = ["x"]
collider_rsq_one = rsquared(collider_H, variables, conditionants).subs( 
        {"beta_z" : 0,
         "beta_y" : 0,
         "beta_x" : 0,
         "beta_yz" : 1,
         "sigma_z" : 1,
         "sigma_x" : 1,
         "sigma_y" : 1})
 

sympy.simplify(collider_rsq_both/collider_rsq_one)


# With confounder

confounder = networkx.DiGraph()
confounder.add_nodes_from(["x", "y", "z"])
confounder.add_edges_from([("x", "y"),
                  ("z", "x"),
                  ("z", "y")])

confounder = populate_attributes(confounder)
confounder_H = calculate_H(confounder)

variables = ["y"]
conditionants = ["x", "z"]
confounder_rsq_both = rsquared(confounder_H, variables, conditionants).subs( 
        {"beta_z" : 0,
         "beta_y" : 0,
         "beta_x" : 0,
         "beta_zy" : 1,
         "beta_zx" : 1,
         "sigma_z" : 1,
         "sigma_x" : 1,
         "sigma_y" : 1})
    

variables = ["y"]
conditionants = ["x"]
confounder_rsq_one = rsquared(confounder_H, variables, conditionants).subs( 
        {"beta_z" : 0,
         "beta_y" : 0,
         "beta_x" : 0,
         "beta_zy" : 1,
         "beta_zx" : 1,
         "sigma_z" : 1,
         "sigma_x" : 1,
         "sigma_y" : 1})    
    

# With nothing

confounder = networkx.DiGraph()
confounder.add_nodes_from(["x", "y", "z"])
confounder.add_edges_from([("x", "y"),
                  ("z", "x")])

confounder = populate_attributes(confounder)
confounder_H = calculate_H(confounder)

variables = ["y"]
conditionants = ["x", "z"]
confounder_rsq_both = rsquared(confounder_H, variables, conditionants).subs( 
        {"beta_z" : 0,
         "beta_y" : 0,
         "beta_x" : 0,
         "beta_zy" : 1,
         "beta_zx" : 1,
         "sigma_z" : 1,
         "sigma_x" : 1,
         "sigma_y" : 1})
    

variables = ["y"]
conditionants = ["x"]
confounder_rsq_one = rsquared(confounder_H, variables, conditionants).subs( 
        {"beta_z" : 0,
         "beta_y" : 0,
         "beta_x" : 0,
         "beta_zx" : 1,
         "sigma_z" : 1,
         "sigma_x" : 1,
         "sigma_y" : 1})      
    
    
# With magic

collider = networkx.DiGraph()
collider.add_nodes_from(["t", "x", "y", "u", "v"])

collider.add_edges_from([("t", "y"),
                  ("v", "x"),
                  ("v", "y"),
                  ("u", "x"),
                  ("u", "t")])

collider = populate_attributes(collider)
collider_H = calculate_H(collider)

variables = ["y"]
conditionants = ["x", "t"]
rsq = rsquared(collider_H, variables, conditionants)
collider_rsq_both = rsquared(collider_H, variables, conditionants).subs( 
        {"beta_z" : 0,
         "beta_y" : 0,
         "beta_x" : 0,
         "beta_yz" : 1,
         "sigma_z" : 1,
         "sigma_x" : 1,
         "sigma_y" : 1})    