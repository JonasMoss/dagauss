# -*- coding: utf-8 -*-
"""
Created on Sat Nov 17 10:36:44 2018

@author: jonas
"""

import networkx
from itertools import product

G = networkx.DiGraph()
G.add_nodes_from(["v", "w", "x", "y", "z"])
G.add_edges_from([("v", "x"),
                  ("w", "y"),
                  ("v", "z"),
                  ("z", "w"),
                  ("z", "x"),
                  ("z", "y"),
                  ("x", "y")])
to_dag(G)

responses = ["x", "y"]
covariates = ["w", "z"]
conditionants = ["v"]
beta(G, responses, covariates + conditionants)

beta(G, responses, covariates )
covariance(G, responses, covariates)

