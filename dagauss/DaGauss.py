"""
This file contains the class defitions for DaGauss.

The directed graph has the following labels reserved:
    [Actually, it only needs conditioned... Which I can delay.]
    Nodes: sigma, beta, conditioned
    Arrows: beta
    
Check for equivalence?

 ** Arrows between two conditioned nodes are redundant. **

    


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
                1.) Removal of arrows: When an arrow is removed, all direct 
                      association is removed.
                2.) Conditioning: Conditioning changes the structure of the 
                      graph, inducing spurious associations.
                3.) Marginalization: Integrate out the variables not there.
              [This is not unique! How do I make it unique?]

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