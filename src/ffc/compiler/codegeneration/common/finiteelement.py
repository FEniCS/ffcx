"Code generation for finite element"

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-01-23 -- 2007-04-19"
__copyright__ = "Copyright (C) 2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# Modified by Kristian Oelgaard 2007

# FFC fem modules
from ffc.fem.finiteelement import *
from ffc.fem.vectorelement import *
from ffc.fem.projection import *
from ffc.fem.dofmap import *

# FFC code generation common modules
from evaluatebasis import *
from evaluatebasisderivatives import *
from utils import *

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
    code["evaluate_basis"] = ["// Not implemented"]
    #code["evaluate_basis"] = evaluate_basis(element, format)

    # Generate code for evaluate_basis_derivatives (FIXME: not implemented)
#    code["evaluate_basis_derivatives"] = ["// Not implemented"]
#    code["evaluate_basis_derivatives"] = evaluate_basis_derivatives(element, format)

    # Generate code for evaluate_dof
    code["evaluate_dof"] = __generate_evaluate_dof(element, format)

    # Generate code for inperpolate_vertex_values
    code["interpolate_vertex_values"] = __generate_interpolate_vertex_values(element, format)

    # Generate code for num_sub_elements
    code["num_sub_elements"] = "%d" % element.num_sub_elements()

    return code

def __generate_evaluate_dof(element, format):
    "Generate code for evaluate_dof"

    # Generate code as a list of lines
    code = []

    # Generate dof map
    dof_map = DofMap(element)
    
    # Check if evaluate_dof is supported
    if dof_map.dof_components() == None:
        code += [format["comment"]("Not implemented (only for Lagrange elements")]
        code += [format["return"](format["floating point"](0.0))]
        return code

    # Get code formats
    block = format["block"]
    separator = format["separator"]
    floating_point = format["floating point"]

    # Get dof coordinates
    cs = dof_map.dof_coordinates()
    X = block(separator.join([block(separator.join([floating_point(x) for x in c])) for c in cs]))

    # Get dof components
    cs = dof_map.dof_components()
    components = block(separator.join(["%d" % c for c in cs]))

    # Compute number of values
    num_values = 1
    for i in range(element.value_rank()):
        num_values *= element.value_dimension(i)

    code += [format["snippet evaluate_dof"](element.cell_dimension()) % \
             (num_values, dof_map.local_dimension(), X, dof_map.local_dimension(), components)]
    
    return code

def __generate_interpolate_vertex_values(element, format):
    "Generate code for interpolate_vertex_values"

    # Check that we have a scalar- or vector-valued element
    if element.value_rank() > 1:
        return format["comment"]("Not implemented (only for scalars or vectors)")

    # Generate code as a list of declarations
    code = []

    # Set vertices (note that we need to use the FIAT reference cells)
    if element.cell_shape() == LINE:
        vertices = [(-1,), (1,)]
    elif element.cell_shape() == TRIANGLE:
        vertices = [(-1, -1), (1, -1), (-1, 1)]
    elif element.cell_shape() == TETRAHEDRON:
        vertices =  [(-1, -1, -1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]

    # Tabulate basis functions at vertices
    table = element.tabulate(0, vertices)

    # Get vector dimension
    if element.value_rank() == 0:
        for i in range(len(vertices)):
            coefficients = table[0][element.cell_dimension()*(0,)][:, i]
            dof_values = [format["dof values"](j) for j in range(len(coefficients))]
            name = format["vertex values"](i)
            value = inner_product(coefficients, dof_values, format)
            code += [(name, value)]
    else:
        for dim in range(element.value_dimension(0)):
            for i in range(len(vertices)):
                coefficients = table[dim][0][element.cell_dimension()*(0,)][:, i]
                dof_values = [format["dof values"](j) for j in range(len(coefficients))]
                name = format["vertex values"](dim*len(vertices) + i)
                value = inner_product(coefficients, dof_values, format)
                code += [(name, value)]

    return code
