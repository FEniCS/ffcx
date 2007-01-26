"Language and format independent code generation"

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-01-23 -- 2007-01-27"
__copyright__ = "Copyright (C) 2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# FFC codegen modules
from finiteelement import *
from dofmap import *
from form import *

def generate_code(name, elements, dof_maps, form, format):
    "Generate code according to given format"

    code = {}

    # Set name
    code["name"] = name

    # Set number of arguments
    code["num_arguments"] = len(elements)

    # Generate code for finite elements
    for i in range(len(elements)):
        code[("finite_element", i)] = generate_finite_element(elements[i], format)

    # Generate code for dof maps
    for i in range(len(dof_maps)):
        code[("dof_map", i)] = generate_dof_map(dof_maps[i], format)

    # Generate code for form
    code["form"] = generate_form(form, format)

    return code
