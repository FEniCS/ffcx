"Code generator"

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-01-23 -- 2007-01-23"
__copyright__ = "Copyright (C) 2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# FFC codegen modules
from finiteelement import *

def generate_code(elements, format):
    "Generate code according to given format"

    code = {}

    # Generate code for finite elements
    for i in range(len(elements)):
        code[("finite_element", i)] = generate_finite_element(elements[i], format)

    return code
