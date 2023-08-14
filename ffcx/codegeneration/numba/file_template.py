# Code generation format strings for UFC (Unified Form-assembly Code)
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2018.

declaration_pre = """
"""

declaration_post = ""

implementation_pre = """
# This code conforms with the UFC specification version {ufcx_version}
# and was automatically generated by FFCx version {ffcx_version}.
# 
# This code was generated with the following options:
#
{options}

import numba
import numpy as np

# ufcx enums
interval = 10
triangle = 20
quadrilateral = 30
tetrahedron = 40
hexahedron = 50
vertex = 60
prism = 70
pyramid = 80

cell = 0
exterior_facet = 1
interior_facet = 2

ufcx_basix_element = 0
ufcx_mixed_element = 1
ufcx_quadrature_element = 2
ufcx_basix_custom_element = 3

"""

implementation_post = """
"""
