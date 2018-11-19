dagauss: Multivariate Gaussian DAGs
==========================

*Note:* This software is in very early development. It has no stable interface.

This module exists to help you investigate multivariate Gaussian DAGs 
symbolically. 

Take the following simple collider:

![image](docs/images/simple_collider.png)

What happens to the estimate of regression coefficient of `x` on `y` when we
condition on `z`? We know it induces a spurious correlation between 
`x` and `y`, but how strong will this spurious correlation be?

Or take the similar problem with a confounder:

![image](docs/images/simple_confounder.png)

What happens to the estimate of regression coefficient of `x` on `y` when we 
do not condition on `z`?

# Examples

To make to output of the examples more appealing, use
```python
import sympy
sympy.init_printing(use_unicode = True)
```

## A confounder
Let's take a look at the confounder example. To use `dagauss` you should 
first make a directed acyclic graph in `networkx`, then make it into a 
`DaGauss` object with `to_dagauss`.

```python
import networkx
import dagauss as dg 

# Make a DAG with confounder z.

G = networkx.DiGraph()
G.add_nodes_from(["x", "y", "z"])
G.add_edges_from([("x", "y"),
                  ("z", "x"),
                  ("z", "y")])
				  
dg.to_dagauss(G)
```

This creates a DAG where each node contains its own variance and mean, and each
edge from "x" to "y" has a `beta_xy` coefficient, giving the conditional 
regression coefficient of `y ~ x`. Visually it's the same graph as the second 
image above.

Let's use this DAG to make compute the conditional mean and covariance of 
`y` and `x` given `z`. We could plausibly be interested in these since we 
wish to measure the causal effect of `x` on `y`, but this is confounded 
by `z` as it now stands.

```python
responses = ["y"]
covariates = ["x"]
conditionants = ["z"]

# beta gives the regression coefficent beta_xy when z is conditioned on.
dg.beta(G, responses, covariates, conditionants)

# The conditional variance of y and x given z:
dg.covariance(G, responses + covariates, conditionants)
```

We can go further and calculate the conditional correlation and conditional 
R squared.

```python
# The conditional correlation of y and x given z:
dg.correlation(G, responses + covariates, conditionants)[0, 1]

# The rsquared of y ~ x | z.
dg.rsquared(G, responses, covariates, conditionants)
```

But we rarely wish to condition on `z`, as this keeps `z` from exerting its 
natural effect on the outcome variables. What we want to do surgery on the DAG 
itself to get to the causal summaries. 
This involves removing all the edges from the parents of the causal variables 
to the causal variables themselves. In this case, the causal variable is `x` 
and its only parent is `z`. Using `sympy`'s `subs` function, we force 
`beta_zx = 0`.

```python
dg.beta(G, responses, covariates).subs({"beta_zx" : 0})
dg.covariance(G, responses + covariates).subs({"beta_zx" : 0})
dg.correlation(G, responses + covariates).subs({"beta_zx" : 0})[0, 1]
dg.rsquared(G, responses, covariates).subs({"beta_zx" : 0})
```

# Development

The interface of this package *will* change. There will be a `DaGauss` class and
the functions used in the exampels above will be methods of the class. There will 
hopefully be more features, such as functions to analyze the causal properties of DAGs, 
maybe a interface to `STAN`, et cetera.

# Installation
You can install it from the command line using `pip`: 

```bash
pip install https://github.com/JonasMoss/dagauss/archive/master.zip
```

# Documentation
No documentation yet.