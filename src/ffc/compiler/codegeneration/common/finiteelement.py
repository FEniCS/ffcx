"Code generation for finite element"

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-01-23 -- 2007-05-18"
__copyright__ = "Copyright (C) 2007 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Kristian Oelgaard 2007

# FIAT modules
from FIAT.shapes import LINE

# FFC fem modules
from ffc.fem.finiteelement import *
from ffc.fem.vectorelement import *
from ffc.fem.projection import *
from ffc.fem.dofmap import *

from ffc.fem.quadratureelement import *

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

    # Disable code generation for unsupported functions of QuadratureElement,
    # (or MixedElements including QuadratureElements)
    if not True in [isinstance(e, QuadratureElement) for e in element.basis_elements()]:
        # Generate code for evaluate_basis
        code["evaluate_basis"] = evaluate_basis(element, format)

        # Generate code for evaluate_basis_derivatives
        code["evaluate_basis_derivatives"] = evaluate_basis_derivatives(element, format)

        # Generate code for inperpolate_vertex_values
        #code["interpolate_vertex_values"] = __generate_interpolate_vertex_values_old(element, format)
        code["interpolate_vertex_values"] = __generate_interpolate_vertex_values(element, format)
    else:
        code["evaluate_basis"] = format["exception"]("evaluate_basis() is not supported for QuadratureElement")
        code["evaluate_basis_derivatives"] = format["exception"]("evaluate_basis_derivatives() is not supported for QuadratureElement")
        code["interpolate_vertex_values"] = format["exception"]("interpolate_vertex_values() is not supported for QuadratureElement")

    # Generate code for evaluate_dof
    code["evaluate_dof"] = __generate_evaluate_dof(element, format)

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
#    if dof_map.dof_coordinates() == None or dof_map.dof_components() == None:
    if None in dof_map.dof_coordinates() or dof_map.dof_components() == None:
        code += [format["exception"]("evaluate_dof not implemented for this type of element")]
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

def __generate_interpolate_vertex_values_old(element, format):
    "Generate code for interpolate_vertex_values"

    # Check that we have a scalar- or vector-valued element
    if element.value_rank() > 1:
        return format["exception"]("interpolate_vertex_values not implemented for this type of element")

    # Generate code as a list of declarations
    code = []

    # Set vertices
    if element.cell_shape() == LINE:
        vertices = [(0.0,), (1.0,)]
    elif element.cell_shape() == TRIANGLE:
        vertices = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]
    elif element.cell_shape() == TETRAHEDRON:
        vertices =  [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0,0), (0,0, 0,0, 1.0)]

    # Tabulate basis functions at vertices
    table = element.tabulate(0, vertices)

    # Get vector dimension
    if element.value_rank() == 0:
        for v in range(len(vertices)):
            coefficients = table[0][element.cell_dimension()*(0,)][:, v]
            dof_values = [format["dof values"](n) for n in range(len(coefficients))]
            name = format["vertex values"](v)
            value = inner_product(coefficients, dof_values, format)
            code += [(name, value)]
    else:
        for dim in range(element.value_dimension(0)):
            for v in range(len(vertices)):
                coefficients = table[dim][0][element.cell_dimension()*(0,)][:, v]
                dof_values = [format["dof values"](n) for n in range(len(coefficients))]
                name = format["vertex values"](dim*len(vertices) + v)
                value = inner_product(coefficients, dof_values, format)
                code += [(name, value)]

    return code

def __generate_interpolate_vertex_values(element, format):
    "Generate code for interpolate_vertex_values"

    # Check that we have a scalar- or vector-valued element
    if element.value_rank() > 1:
        return format["exception"]("interpolate_vertex_values not implemented for this type of element")

    # Generate code as a list of declarations
    code = []

    # Set vertices (note that we need to use the FIAT reference cells)
    if element.cell_shape() == LINE:
        vertices = [(0,), (1,)]
    elif element.cell_shape() == TRIANGLE:
        vertices = [(0, 0), (1, 0), (0, 1)]
    elif element.cell_shape() == TETRAHEDRON:
        vertices =  [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]

    # Extract nested sub elements
    sub_elements = element.basis_elements()

    # Iterate over sub elements
    offset_dof_values = 0
    offset_vertex_values = 0
    need_jacobian = False

    for sub_element in sub_elements:

        # Tabulate basis functions at vertices
        table = sub_element.tabulate(0, vertices)

        # Check which transform we should use to map the basis functions
        mapping = pick_first([sub_element.value_mapping(dim) for dim in range(sub_element.value_dimension(0))])

        # Generate different code depending on mapping
        if mapping == Mapping.AFFINE:

            code += [format["comment"]("Evaluate at vertices and use affine mapping")]

            # Handle scalars and vectors
            if sub_element.value_rank() == 0:
                for v in range(len(vertices)):
                    coefficients = table[0][sub_element.cell_dimension()*(0,)][:, v]
                    dof_values = [format["dof values"](offset_dof_values + n) for n in range(len(coefficients))]
                    name = format["vertex values"](offset_vertex_values + v)
                    value = inner_product(coefficients, dof_values, format)
                    code += [(name, value)]
            else:
                for dim in range(sub_element.value_dimension(0)):
                    for v in range(len(vertices)):
                        coefficients = table[dim][0][sub_element.cell_dimension()*(0,)][:, v]
                        dof_values = [format["dof values"](offset_dof_values + n) for n in range(len(coefficients))]
                        name = format["vertex values"](offset_vertex_values + dim*len(vertices) + v)
                        value = inner_product(coefficients, dof_values, format)
                        code += [(name, value)]

        elif (mapping == Mapping.CONTRAVARIANT_PIOLA or mapping == Mapping.COVARIANT_PIOLA):

            code += [format["comment"]("Evaluate at vertices and use Piola mapping")]

            # Remember to add code later for Jacobian
            need_jacobian = True

            # Check that dimension matches for Piola transform
            if not sub_element.value_dimension(0) == sub_element.cell_dimension():
                raise RuntimeError, "Vector dimension of basis function does not match for Piola transform."

            # Get entities for the dofs
            dof_entities = DofMap(sub_element).dof_entities()

            for dim in range(sub_element.value_dimension(0)):
                for v in range(len(vertices)):
                    terms = []
                    for n in range(sub_element.space_dimension()):
                        # Get basis function values at vertices
                        coefficients = [table[j][0][sub_element.cell_dimension()*(0,)][n, v] for j in range(sub_element.value_dimension(0))]

                        if mapping == Mapping.COVARIANT_PIOLA:
                            # Get row of inverse transpose Jacobian
                            jacobian_row = [format["transform"](Transform.JINV, j, dim, None) for j in range(sub_element.cell_dimension())]
                        else:
                            # mapping == Mapping.CONTRAVARIANT_PIOLA:
                            # Get row of Jacobian
                            jacobian_row = [format["transform"](Transform.J, j, dim, None) for j in range(sub_element.cell_dimension())]
                            
                        # Multiply vector-valued basis function with Jacobian
                        basis_function = inner_product(coefficients, jacobian_row, format)
                        # Add paranthesis if necessary
                        if "+" in basis_function or "-" in basis_function: # Cheating, should use dictionary
                            basis_function = format["grouping"](basis_function)
                            # Multiply with dof value
                        factors = [format["dof values"](offset_dof_values + n), basis_function]
                        # Add term
                        if not basis_function == format["floating point"](0):
                            terms += [format["multiply"](factors)]
                    if len(terms) > 1:
                        sum = format["grouping"](format["add"](terms))
                    else:
                        sum = format["add"](terms)
                    name = format["vertex values"](offset_vertex_values + dim*len(vertices) + v)
                    if mapping == Mapping.CONTRAVARIANT_PIOLA:
                        value = format["multiply"]([format["inverse"](format["determinant"](None)), sum])
                    else: 
                        value = format["multiply"]([sum])
                    code += [(name, value)]

        else:
            raise RuntimeError, "Unknown mapping: " + str(mapping)

        offset_dof_values    += sub_element.space_dimension()
        offset_vertex_values += len(vertices)*sub_element.value_dimension(0)

    # Insert code for computing quantities needed for Piola mapping
    if need_jacobian:
        code.insert(0, format["snippet jacobian"](element.cell_dimension()) % {"restriction": ""})        
    
    return code
