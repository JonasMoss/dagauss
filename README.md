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
import networkx as nx

# Add a DAG with confounder z.

G = nx.DiGraph()
G.add_nodes_from(["x", "y", "z"])
G.add_edges_from([("x", "y"),
                  ("z", "x"),
                  ("z", "y")])

# Calculate its non-DAG:
G = populate_attributes(G)
H = calculate_H(G)

# Find its variance and covariance.
mean(H)
covariance(H)
```

# Installation
Copy the scripts and use them.

# Documentation
No documentation yet.