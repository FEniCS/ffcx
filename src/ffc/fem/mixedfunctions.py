"Factory functions for mixed basis functions and functions"

_author__ = "Anders Logg (logg@simula.no)"
__date__ = "2005-09-16 -- 2007-03-29"
__copyright__ = "Copyright (C) 2005-2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# Modified by Marie E. Rognes (meg@math.uio.no) 2007

# FFC compiler.language modules
from ffc.compiler.language.algebra import *

# FFC fem modules
from ffc.fem.mixedelement import *

def BasisFunctions(element, functiontype = BasisFunction):
    "Create tuple of BasisFunctions from given MixedElement"
    if not isinstance(element, MixedElement):
        raise RuntimeError, "Basis function tuple must be created from mixed element."
    # Create basis functions for mixed element
    basisfunctions = tabulate_components(element, functiontype)
    return basisfunctions


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
    # Create function for mixed element
    functions = tabulate_components(element, Function)
    return functions

def tabulate_components(element, functiontype):
    "Explicitely tabulate components of mixed (basis) functions"
    vector = functiontype(element)
    components = []
    offset = 0
    # Pick components
    for i in range(element.num_sub_elements()):
        sub_element = element.sub_element(i)
        if sub_element.value_rank() == 0:
            sub_vector = vector[Index(offset)]
            offset += 1
        elif sub_element.value_rank() == 1:
            sub_vector = [vector[Index(k)] for k in range(offset, sub_element.value_dimension(0)+offset)]
            offset += sub_element.value_dimension(0)
        else:
            raise RuntimeError, "Mixed elements can only be created from scalar or vector-valued elements."
        components += [sub_vector]
    return tuple(components)
