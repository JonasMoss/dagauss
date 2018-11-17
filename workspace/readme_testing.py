# -*- coding: utf-8 -*-
"""
Created on Sat Nov 17 20:58:22 2018

@author: jonas
"""

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

# The rsquared of y ~ x | z.
rsquared(G, responses, covariates, conditionants)