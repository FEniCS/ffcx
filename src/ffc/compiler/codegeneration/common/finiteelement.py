"Code generation for finite element"

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-01-23 -- 2007-02-06"
__copyright__ = "Copyright (C) 2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# FFC fem modules
from ffc.fem.finiteelement import *

# FFC evaluatebasis module
from ffc.compiler.codegeneration.common.evaluatebasis import *

def generate_finite_element(element, format):
    """Generate dictionary of code for the given finite element
    according to the given format"""

    code = {}

    # Generate code for signature
    code["signature"] = element.signature()

    # Generate code for cell_shape
    code["cell_shape"] = format["cell shape"](element.cell_shape())
    
    # Generate code for space_dimension
    code["space_dimension"] = "%d" % element.space_dimension()

    # Generate code for value_rank
    code["value_rank"] = "%d" % element.value_rank()

    # Generate code for value_dimension
    code["value_dimension"] = ["%d" % element.value_dimension(i) for i in range(max(element.value_rank(), 1))]

    # Generate code for evaluate_basis (FIXME: not implemented)
#    code["evaluate_basis"] = "// Not implemented"
    code["evaluate_basis"] = evaluate_basis(element, format)

    # Generate code for evaluate_dof (FIXME: not implemented)
    code["evaluate_dof"] = "// Not implemented\nreturn 0.0;"

    # Generate code for inperpolate_vertex_values (FIXME: not implemented)
    code["interpolate_vertex_values"] = "// Not implemented"

    # Generate code for num_sub_elements
    code["num_sub_elements"] = "%d" % element.num_sub_elements()

    return code
