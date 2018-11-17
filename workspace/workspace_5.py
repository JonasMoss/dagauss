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