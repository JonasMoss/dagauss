# -*- coding: utf-8 -*-
"""
Created on Sat Nov 17 20:58:22 2018

@author: jonas
"""

G = networkx.DiGraph()
G.add_nodes_from(["x", "y"])
G.add_edges_from([("x", "y")])
to_dag(G)

responses = ["y"]
covariates = ["x"]
conditionants = []

subs =  {"beta_x"  : 0,
         "beta_y"  : 0,
         "beta_xy" : 1,
         "sigma_x" : 1,
         "sigma_y" : 1}

rsqsign(G, responses, covariates)
correlation(G, responses, covariates)