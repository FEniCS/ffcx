"Factory functions for mixed basis functions and functions"

_author__ = "Anders Logg (logg@simula.no)"
__date__ = "2005-09-16 -- 2007-03-29"
__copyright__ = "Copyright (C) 2005-2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# FFC compiler.language modules
from ffc.compiler.language.algebra import *

# FFC fem modules
from ffc.fem.mixedelement import *

def BasisFunctions(element, functiontype = BasisFunction):
    "Create tuple of BasisFunctions from given MixedElement"
    if not isinstance(element, MixedElement):
        raise RuntimeError, "Basis function tuple must be created from mixed element."
    # Create basis function for mixed element
    vector = functiontype(element)
    # Pick components/subvectors of the mixed basis function
    sub_vectors = []
    offset = 0
    for i in range(element.num_sub_elements()):
        sub_element = element.sub_element(i)
        if sub_element.value_rank() == 0:
            sub_vector = vector.pick_component_default(offset)
            offset += 1
        elif sub_element.value_rank() == 1:
            if sub_element.mapping() == "Piola":
                sub_vector = [vector.pick_component_piola(k) for k in range(0, sub_element.value_dimension(0))]
            else:
                sub_vector = [vector.pick_component_default(i) for i in range(offset, offset + sub_element.value_dimension(0))]
            offset += sub_element.value_dimension(0)
        else:
            raise RuntimeError, "Mixed elements can only be created from scalar or vector-valued elements."
        sub_vectors += [sub_vector]
    return tuple(sub_vectors)

def TestFunctions(element):
    "Create tuple of TestFunctions from given MixedElement"
    return BasisFunctions(element, TestFunction)

def TrialFunctions(element):
    "Create tuple of TrialFunctions from given MixedElement"
    return BasisFunctions(element, TrialFunction)

def Functions(element):
    "Create tuple of Functions from given MixedElement."
    if not isinstance(element, MixedElement):
        raise RuntimeError, "Function tuple must be created from mixed element."
    # Create function fox mixed element
    vector = Function(element)
    # Pick components/subvectors of the mixed basis function
    subvectors = []
    offset = 0
    for e in element.elements:
        if e.value_rank() == 0:
            subvector = vector[offset]
            offset += 1
        elif e.value_rank() == 1:
            subvector = [vector[i] for i in range(offset, offset + e.value_dimension(0))]
            offset += e.value_dimension(0)
        else:
            raise RuntimeError, "Mixed elements can only be created from scalar or vector-valued elements."
        subvectors += [subvector]
    return tuple(subvectors)
