"Code generation for finite element"

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-01-23 -- 2007-01-23"
__copyright__ = "Copyright (C) 2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# FFC fem modules
from ffc.fem.finiteelement import *

def generate_finite_element(element, format):
    """Generate dictionary of code for the given finite element
    according to the given format."""

    code = {}

    # Generate code for signature
    code["signature"] = element.signature()

    # Generate code for cell_shape
    code["cell_shape"] = "return ufc::%s;" % shape_to_string[element.cell_shape()]
    
    # Generate code for space_dimension
    code["space_dimension"] = "%d" % element.space_dimension()

    # Generate code for value_rank
    code["value_rank"] = "%d" % element.value_rank()

    # Generate code for value_dimension
    if element.value_rank() == 0:
        body = "return %d;" % element.value_dimension(0)
    else:
        body = "switch ( i )\n{\n"
        for i in range(element.value_rank()):
            body += "case %d:\n  return %d;\n  break;\n" % element.value_dimension(i)
        body += "default:\n  return 0;\n}"
    code["value_dimension"] = body

    # Generate code for evaluate_basis (FIXME: not implemented)
    code["evaluate_basis"] = "// Not implemented"

    # Generate code for evaluate_dof (FIXME: not implemented)
    code["evaluate_dof"] = "// Not implemented\nreturn 0.0;"

    # Generate code for inperpolate_vertex_values (FIXME: not implemented)
    code["interpolate_vertex_values"] = "// Not implemented"

    # Generate code for num_sub_elements
    code["num_sub_elements"] = "return %d;" % element.num_sub_elements()

    return code
