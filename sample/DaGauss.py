"""
This file contains the class defitions for DaGauss.

The directed graph has the following labels reserved:
    [Actually, it only needs conditioned... Which I can delay.]
    Nodes: sigma, beta, conditioned
    Arrows: beta
    


Prototype:
    * A directed acyclic graph G. If a node has the attribute 
        conditioned = True, it is conditioned on.
    * Graph attribute: A non-directed graph with the dependency structure.
    * Graph attribute: The mean vector and covariance matrix.
    * Methods: 
        Calculate conditional means and covariances.

        Transform. Using a subgraph of self as input, find the induced 
                * Check that H is a subgraph of self.
            DAG. This is done in two steps:
                1.) Removal of arrows. When an arrow is removed, all direct 
                    association is removed.
                2.) Conditioning. 
                3.) Marginalization.
            done by marginalization of nodes left out and 

"""

from itertools import product
import sympy
import networkx 

class DaGauss(networkx.DiGraph):

    def __init__(self):
        self.firstname = first
        self.lastname = last

    def Name(self):
        return self.firstname + " " + self.lastname