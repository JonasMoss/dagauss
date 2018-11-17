# -*- coding: utf-8 -*-
"""
Created on Sat Nov 17 23:01:41 2018

@author: jonas
"""

N = sympy.Matrix([[0, 1], [0, 0]])



M = sympy.Matrix([[1, 1], 
                  [ 1, 1]])
Y = M

(I + Y).eigenvals()
((I + Y)*Y*(I + Y).T - Y).eigenvals()
