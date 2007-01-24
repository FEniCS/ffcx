"Code generator"

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-01-23 -- 2007-01-23"
__copyright__ = "Copyright (C) 2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# FFC codegen modules
from finiteelement import *
from dofmap import *

def generate_code(elements, dofmaps, format):
    "Generate code according to given format"

    code = {}

    # Generate code for finite elements
    for i in range(len(elements)):
        code[("finite_element", i)] = generate_finite_element(elements[i], format)

    # Generate code for dof maps
    for i in range(len(dofmaps)):
        code[("dof_map", i)] = generate_dof_map(dofmaps[i], format)

    return code
