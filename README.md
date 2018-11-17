dagauss: Multivariate Gaussian DAGs
==========================

*Note:* Not even a pre-release, this is. 

Use this module to investigate multivariate Gaussian DAGs symbolically. Take the following simple collider:

![image](docs/images/simple_collider.png)

What happens to the estimate of regression coefficient of `x` on `y` when we condition on `z`?

Or take the similar problem with a confounder:

![image](docs/images/simple_confounder.png)

What happens to the estimate of regression coefficient of `x` on `y` when we do not condition on `z`?

# Example

```python
import networkx

# Add a DAG with confounder z.

G = networkx.DiGraph()
G.add_nodes_from(["x", "y", "z"])
G.add_edges_from([("x", "y"),
                  ("z", "x"),
                  ("z", "y")])
to_dag(G)

# Find some causal summaries of x on beta.

responses = ["y"]
covariates = ["x"]
conditionants = ["z"]

# beta gives the regression coefficent beta_xy when z is conditioned on.
beta(G, responses, covariates, conditionants)

# The conditional variance of y and x given z:
covariance(G, responses + covariates, conditionants)

# The R^2 of y ~ x | z.
rsquared(G, responses, covariates, conditionants)
```

# Installation
Copy the scripts and use them.

# Documentation
No documentation yet.