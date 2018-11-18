    

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