_author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-03-20 -- 2007-03-20"
__copyright__ = "Copyright (C) 2007 Anders Logg"
__license__  = "GNU GPL Version 2"

from finiteelement import *
from mixedelement import *

def VectorElement(family, shape, degree, num_sub_elements=None):
    "Create vector finite element as the repeated tensor product of a given element"

    # Choose default number of sub elements
    num_sub_elements = shape_to_dim[string_to_shape[shape]]

    # Create mixed element by repeating the "scalar" element
    return MixedElement([FiniteElement(family, shape, degree) for i in range(num_sub_elements)])
