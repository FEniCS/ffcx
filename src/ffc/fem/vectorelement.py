"Factory function for generating vector elements (and vector QuadratureElements) as mixed elements"

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-03-20 -- 2008-01-16"
__copyright__ = "Copyright (C) 2007 - 2008 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Kristian Oelgaard 2008

from finiteelement import *
from quadratureelement import *
from mixedelement import *

def VectorElement(family, shape, degree, num_sub_elements=None):
    "Create vector finite element as the repeated tensor product of a given element"

    # Choose default number of sub elements if number is not specified by user
    if not num_sub_elements:
        num_sub_elements = shape_to_dim[string_to_shape[shape]]

    # Create mixed element by repeating the "scalar" element
    return MixedElement([FiniteElement(family, shape, degree) for i in range(num_sub_elements)])

def VectorQuadratureElement(shape, degree, num_sub_elements=None):
    "Create vector quadrature element as the repeated tensor product of a given element"

    # Choose default number of sub elements if number is not specified by user
    if not num_sub_elements:
        num_sub_elements = shape_to_dim[string_to_shape[shape]]

    # Create mixed element by repeating the "scalar" element
    return MixedElement([QuadratureElement(shape, degree) for i in range(num_sub_elements)])

