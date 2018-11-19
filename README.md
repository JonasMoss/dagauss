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
import dagauss as dg 

# Make a DAG with confounder z.

G = networkx.DiGraph()
G.add_nodes_from(["x", "y", "z"])
G.add_edges_from([("x", "y"),
                  ("z", "x"),
                  ("z", "y")])

# Make G into a dagauss object:
dg.to_dagauss(G)

# Find some causal summaries of x on beta.

responses = ["y"]
covariates = ["x"]
conditionants = ["z"]

# beta gives the regression coefficent beta_xy when z is conditioned on.
dg.beta(G, responses, covariates, conditionants)

# The conditional variance of y and x given z:
dg.covariance(G, responses + covariates, conditionants)

# The rsquared of y ~ x | z.
dg.rsquared(G, responses, covariates, conditionants)
```

The interface of this package will change. There will be a `DaGauss` class and
the functions above will be methods of the class.

# Installation
You can install it from the command line.

```bash
pip install https://github.com/django/django/archive/master.zip
```

# Documentation
No documentation yet.