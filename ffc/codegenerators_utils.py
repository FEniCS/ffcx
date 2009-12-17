"""This module contains all the helper functions for the codegenerators module
and the functions from the old utils and codeutils modules."""

__author__ = "Anders Logg (logg@simula.no) and Kristian B. Oelgaard (k.b.oelgard@gmail.com)"
__date__ = "2009-12-09"
__copyright__ = "Copyright (C) 2009 Anders Logg and Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# Last changed: 2009-12-17

# Python modules.
import numpy

# UFL modules.
from ufl import FiniteElement as UFLFiniteElement

# FFC modules.
from utils import pick_first
from log import error

#------------------------------------------------------------------------------
# From utils.py.
def indent(s, n):
    "Indent each row of the given string s with n spaces"
    indentation = " "*n
    return indentation + ("\n" + indentation).join(s.split("\n"))
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# From codeutils.py.
def inner_product(a, b, format):
    """Generate code for inner product of a and b, where a is a list
    of floating point numbers and b is a list of symbols."""

    # Check input
    if not len(a) == len(b):
        error("Dimensions don't match for inner product.")

    # Prefetch formats to speed up code generation
    format_add            = format["add"]
    format_subtract       = format["subtract"]
    format_multiply       = format["multiply"]
    format_floating_point = format["floating point"]
    format_epsilon        = format["epsilon"]

    # Add all entries
    value = None
    for i in range(len(a)):

        # Skip terms where a is almost zero
        if abs(a[i]) <= format_epsilon:
            continue

        # Fancy handling of +, -, +1, -1
        if value:
            if abs(a[i] - 1.0) < format_epsilon:
                value = format_add([value, b[i]])
            elif abs(a[i] + 1.0) < format_epsilon:
                value = format_subtract([value, b[i]])
            elif a[i] > 0.0:
                value = format_add([value, format_multiply([format_floating_point(a[i]), b[i]])])
            else:
                value = format_subtract([value, format_multiply([format_floating_point(-a[i]), b[i]])])
        else:
            if abs(a[i] - 1.0) < format_epsilon or abs(a[i] + 1.0) < format_epsilon:
                value = b[i]
            else:
                value = format_multiply([format_floating_point(a[i]), b[i]])

    return value or format_floating_point(0.0)

def tabulate_matrix(matrix, format):
    "Function that tabulates the values of a matrix, into a two dimensional array."

    # Check input
    if not len(numpy.shape(matrix)) == 2:
        error("This is not a matrix.")

    # Prefetch formats to speed up code generation
    format_block          = format["block"]
    format_separator      = format["separator"]
    format_floating_point = format["floating point"]
    format_epsilon        = format["epsilon"]

    # Get size of matrix
    num_rows = numpy.shape(matrix)[0]
    num_cols = numpy.shape(matrix)[1]

    # Set matrix entries equal to zero if their absolute values is smaller than format_epsilon
    for i in range(num_rows):
        for j in range(num_cols):
            if abs(matrix[i][j]) < format_epsilon:
                matrix[i][j] = 0.0

    # Generate array of values
    value = format["new line"] + format["block begin"]
    rows = []

    for i in range(num_rows):
        rows += [format_block(format_separator.join([format_floating_point(matrix[i,j])\
                 for j in range(num_cols)]))]

    value += format["block separator"].join(rows)
    value += format["block end"]

    return value

class IndentControl:
    "Class to control the indentation of code"

    def __init__(self):
        "Constructor"
        self.size = 0
        self.increment = 2

    def increase(self):
        "Increase indentation by increment"
        self.size += self.increment

    def decrease(self):
        "Decrease indentation by increment"
        self.size -= self.increment

    def indent(self, a):
        "Indent string input string by size"
        return indent(a, self.size)
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# From finiteelement.py
def __extract_sub_elements(element, parent):
    """Recursively extract sub elements as a list of tuples where
    each tuple consists of a tuple labeling the sub element and
    the sub element itself."""

    if element.num_sub_elements() == 1:
        return [(parent, element)]
    sub_elements = []
    for i in range(element.num_sub_elements()):
        sub_elements += __extract_sub_elements(element.sub_element(i), parent + (i,))
    return sub_elements + [(parent, element)]

# FIXME: This is C++ dependent, move relevant parts to ufc_format.py
def __generate_evaluate_dof(element, format):
    "Generate code for evaluate_dof"

    # Generate code as a list of lines
    code = []

    # Get code formats
    block = format["block"]
    separator = format["separator"]
    floating_point = format["floating point"]
    comment = format["comment"]

    # Generate dof map and get _copy_ of dof representations.
    dof_map = DofMap(element)
    dofs = [DofRepresentation(dof) for dof in dof_map.dual_basis()]
    num_dofs = len(dofs)

    # For ease in the code generalization, pad the points, directions
    # and weights with zeros according to the maximal number of
    # points. (Hence, the copy of the dofs above.)
    max_num_points = dof_map.get_max_num_of_points()
    num_points_per_dof = [dof.pad_points_and_weights(max_num_points)
                          for dof in dofs]

    # Compute the value dimension of the functions
    num_values = 1
    for i in range(element.value_rank()):
        num_values *= element.value_dimension(i)
    # Check that the value dimension is the same for all dofs and that
    # it matches the dimension of the function
    value_dim = pick_first([dof.value_dim() for dof in dofs])
    if value_dim != num_values:
        error("Directional component does not match vector dimension")

    # Initialize the points, weights and directions for the dofs:
    code += [comment("The reference points, direction and weights:")]
    s = block(separator.join(
        [block(separator.join([block(separator.join([floating_point(c)
                                                     for c in point]))
                               for point in dof.points]))
         for dof in dofs]))
    code  += ["%sX[%d][%d][%d] = %s;" % (format["table declaration"],
                                         num_dofs, max_num_points,
                                         element.geometric_dimension(), s)]
    s = block(separator.join(
        [block(separator.join([floating_point(w) for w in dof.weights]))
         for dof in dofs]))
    code += ["%sW[%d][%d] = %s;" % (format["table declaration"],
                                    num_dofs, max_num_points, s)]

    s = block(separator.join(
        [block(separator.join([block(separator.join([floating_point(dk)
                                                    for dk in d]))
                               for d in dof.directions]))
         for dof in dofs]))
    code += ["%sD[%d][%d][%d] = %s;" % (format["table declaration"],
                                        num_dofs, max_num_points,
                                        num_values, s)]
    code += [""]

    # Compute the declarations needed for function mapping and the
    # code for mapping each set of function values:
    (map_declarations, map_values_code) = \
                       __map_function_values(num_values, element, format)
    code += [map_declarations]

    # Loop over the number of points (if more than one) and evaluate
    # the functional
    code += ["%sresult = 0.0;" % format["float declaration"]]
    code += [comment("Iterate over the points:") ]
    tab = 0
    endloop = ""

    # If there is more than one point, we need to add a table of the
    # number of points per dof, add the loop, add indentation and add
    # an end brackets
    index = "0"
    if max_num_points > 1:
        num_points_per_dof_code = block(separator.join([str(n) for n in num_points_per_dof]))
        code += ["%sns[%d] = %s;" % (format["static const uint declaration"],
                                      num_dofs, num_points_per_dof_code)]
        code += ["%s%s" % (format["loop"]("j", "0", "ns[i]"), " {")]
        (tab, endloop, index) = (2, "\n} // End for", "j")

    # Map the points from the reference onto the physical element
    code += [indent(format["snippet map_onto_physical"](element.geometric_dimension())
                    % {"j": index}, tab)]

    # Evaluate the function at the physical points
    code += [indent(comment("Evaluate function at physical points"), tab)]
    code += [indent("double values[%d];" % num_values, tab)]
    code += [indent("f.evaluate(values, y, c);\n", tab)]

    # Map the function values according to the given mapping(s)
    code += [indent(comment("Map function values using appropriate mapping"),
                    tab)]
    code += [indent(map_values_code, tab)]
    code += [""]

    # Note that we do not map the weights (yet).
    code += [indent(comment("Note that we do not map the weights (yet)."),tab)]
    code += [""]

    # Take the directional components of the function values and
    # multiply by the weights:
    code += [indent(format["snippet calculate dof"] % {"dim": value_dim,
                                                       "index": index}, tab)]
    # End possible loop
    code += [endloop]

    # Return the calculated value
    code += [format["return"]("result")]
    return code

# FIXME: This is C++ dependent, move relevant parts to ufc_format.py
def __map_function_values(num_values, element, format):

    block = format["block"]
    separator = format["separator"]
    comment = format["comment"]
    precode = []
    code = []

    # Compute the mapping associated with each dof:
    # Looping the extracted element should have the same effect as using the
    # old element.space_mapping() function.
    # mappings = [mapping_to_int[element.space_mapping(i)]
    #             for i in range(element.space_dimension())]
    mappings = []
    for e in element.extract_elements():
        for d in range(e.space_dimension()):
            mappings.append(mapping_to_int[e.mapping()])

    whichmappings = set(mappings)

    # If there is more than one mapping involved, we will need to
    # keep track of them at runtime:
    if len(whichmappings) > 1:
        precode += ["%smappings[%d] = %s;" %
                    (format["static const uint declaration"], len(mappings),
                     block(separator.join(([str(m) for m in mappings]))))]

    # Check whether we will need a piola
    contrapiola_present = mapping_to_int[CONTRAVARIANT_PIOLA] in whichmappings
    copiola_present = mapping_to_int[COVARIANT_PIOLA] in whichmappings
    piola_present = contrapiola_present or copiola_present

    # Add code for the jacobian if we need it for the
    # mappings. Otherwise, just add code for the vertex coordinates
    if contrapiola_present:
        # If contravariant piola: Will need J, det J and J^{-1}
        precode += [format["snippet jacobian"](element.geometric_dimension())
                 % {"restriction":""}]
        precode += ["\ndouble copyofvalues[%d];" % num_values]
    elif copiola_present:
        # If covariant piola: Will need J only
        precode += [format["snippet only jacobian"](element.geometric_dimension())
                 % {"restriction":""}]
        precode += ["\ndouble copyofvalues[%d];" % num_values]
    else:
        precode += [format["get cell vertices"]]

    # We have to add offsets to the code if there are mixed
    # piola-mapped elements with an offset. (Ex: DG0 + RT)
    offset = ""
    if element.num_sub_elements() > 1 and piola_present:
        value_offsets = []
        adjustment = 0
        for i in range(element.num_sub_elements()):
            subelement = element.sub_element(i)
            value_offsets += [adjustment]*subelement.space_dimension()
            adjustment += subelement.value_dimension(0)

        # if mapping[i] != AFFINE == 0 and value_offset[i] !=
        # 0 for all i, then need_offsets = True:
        need_offsets = bool(max([mappings[i] and value_offsets[i]
                                 for i in range(len(mappings))]))
        if need_offsets:
            precode += ["const int offsets[%d] = %s;" %
                        (len(value_offsets),
                         block(separator.join([str(o) for o in value_offsets])))]
            offset = "offsets[i] + "


    # Then it just remains to actually add the different mappings to the code:
    n = element.geometric_dimension()
    mappings_code = {mapping_to_int[AFFINE]: __affine_map(),
                     mapping_to_int[CONTRAVARIANT_PIOLA]:
                     __contravariant_piola(n, offset),
                     mapping_to_int[COVARIANT_PIOLA]: __covariant_piola(n, offset)}

    ifs = ["mappings[i] == %d" % mapping for mapping in whichmappings]
    cases = [mappings_code[mapping] for mapping in whichmappings]
    code  += [__generate_if_block(ifs, cases,
                                  comment("Other mappings not applicable."))]
    return ("\n".join(precode), "\n".join(code))

# FIXME: This is C++ dependent, move to ufc_format.py
def __affine_map():
    return "// Affine map: Do nothing"

# FIXME: This is C++ dependent, move to ufc_format.py
def __contravariant_piola(dim, offset=""):
    code = []
    code += ["// Copy old values:"]
    for i in range(dim):
        code += ["copyofvalues[%s%d] = values[%s%d];" % (offset, i, offset, i)]
    code += ["// Do the inverse of div piola "]
    for i in range(dim):
        terms = ["Jinv_%d%d*copyofvalues[%s%d]" % (i, j, offset, j)
                 for j in range(dim)]
        code += ["values[%s%d] = detJ*(%s);" % (offset, i, "+".join(terms))]
    return "\n".join(code)

# FIXME: This is C++ dependent, move to ufc_format.py
def __covariant_piola(dim, offset=""):
    code = []
    code += ["// Copy old values:"]
    for i in range(dim):
        code += ["copyofvalues[%s%d] = values[%s%d];" % (offset, i, offset, i)]
    code += ["// Do the inverse of curl piola "]
    for i in range(dim):
        terms = ["J_%d%d*copyofvalues[%s%d]" % (j, i, offset, j)
                 for j in range(dim)]
        code += ["values[%s%d] = %s;" % (offset, i, "+".join(terms))]
    return  "\n".join(code)

# FIXME: This is C++ dependent, move to ufc_format.py
def __generate_if_block(ifs, cases, default = ""):
    "Generate if block from given ifs and cases"
    if len(ifs) != len(cases):
        error("Mismatch of dimensions ifs and cases")

    # Special case: no cases
    if len(ifs) == 0:
        return default

    # Special case: one case
    if len(ifs) == 1:
        return cases[0]

    # Create ifs
    code = "if (%s) { \n" % ifs[0]
    code += "%s\n" % indent(cases[0], 2)
    for i in range(len(ifs)-1):
        code += "} else if (%s) {\n %s \n" % (ifs[i+1], indent(cases[i+1], 2))
    code += "} else { \n %s \n}" % indent(default,2)
    return code

def __generate_interpolate_vertex_values(element, format):
    "Generate code for interpolate_vertex_values"

    # Check that we have a scalar- or vector-valued element
    if element.value_rank() > 1:
        return format["exception"]("interpolate_vertex_values not implemented for this type of element")

    # Generate code as a list of declarations
    code = []

    # Get reference cell vertices
    vertices = get_vertex_coordinates(element.cell_domain())

#    # Set vertices (note that we need to use the FIAT reference cells)
#    if element.cell().domain() == "interval":
#        vertices = [(0,), (1,)]
#    elif element.cell().domain() == "triangle":
#        vertices = [(0, 0), (1, 0), (0, 1)]
#    elif element.cell().domain() == "tetrahedron":
#        vertices =  [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]

    # Extract nested sub elements
    sub_elements = element.extract_elements()

    # Iterate over sub elements
    offset_dof_values = 0
    offset_vertex_values = 0
    need_jacobian = False

    # Compute data size
    size = 0
    for sub_element in sub_elements:
        size += sub_element.value_dimension(0)

    for sub_element in sub_elements:

        # Tabulate basis functions at vertices
        table = sub_element.tabulate(0, vertices)

        # Check which transform we should use to map the basis functions
        mapping = sub_element.mapping()

        # Generate different code depending on mapping
        if mapping == AFFINE:

            code += [format["comment"]("Evaluate at vertices and use affine mapping")]

            # Handle scalars and vectors
            if sub_element.value_rank() == 0:
                for v in range(len(vertices)):
                    coefficients = table[0][sub_element.cell().geometric_dimension()*(0,)][:, v]
                    dof_values = [format["dof values"](offset_dof_values + n) for n in range(len(coefficients))]
                    name = format["vertex values"](size*v + offset_vertex_values)
                    value = inner_product(coefficients, dof_values, format)
                    code += [(name, value)]
            else:
                for dim in range(sub_element.value_dimension(0)):
                    for v in range(len(vertices)):
                        coefficients = table[dim][0][sub_element.cell().geometric_dimension()*(0,)][:, v]
                        dof_values = [format["dof values"](offset_dof_values + n) for n in range(len(coefficients))]
                        name = format["vertex values"](size*v + offset_vertex_values + dim)
                        value = inner_product(coefficients, dof_values, format)
                        code += [(name, value)]

        elif (mapping == CONTRAVARIANT_PIOLA or mapping == COVARIANT_PIOLA):

            code += [format["comment"]("Evaluate at vertices and use Piola mapping")]

            # Remember to add code later for Jacobian
            need_jacobian = True

            # Check that dimension matches for Piola transform
            if not sub_element.value_dimension(0) == sub_element.geometric_dimension():
                error("Vector dimension of basis function does not match for Piola transform.")

            # Get entities for the dofs
            dof_entities = DofMap(sub_element).dof_entities()

            for dim in range(sub_element.value_dimension(0)):
                for v in range(len(vertices)):
                    terms = []
                    for n in range(sub_element.space_dimension()):
                        # Get basis function values at vertices
                        coefficients = [table[j][0][sub_element.cell().geometric_dimension()*(0,)][n, v] for j in range(sub_element.value_dimension(0))]

                        if mapping == COVARIANT_PIOLA:
                            # Get row of inverse transpose Jacobian
                            jacobian_row = [format["transform"]("JINV", j, dim, None) for j in range(sub_element.cell().geometric_dimension())]
                        else:
                            # mapping == CONTRAVARIANT_PIOLA:
                            # Get row of Jacobian
                            jacobian_row = [format["transform"]("J", j, dim, None) for j in range(sub_element.cell().geometric_dimension())]

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
                    name = format["vertex values"](size*v + offset_vertex_values + dim)
                    if mapping == CONTRAVARIANT_PIOLA:
                        value = format["multiply"]([format["inverse"](format["determinant"](None)), sum])
                    else:
                        value = format["multiply"]([sum])
                    code += [(name, value)]

        else:
            error("Unknown mapping: " + str(mapping))

        offset_dof_values    += sub_element.space_dimension()
        offset_vertex_values += sub_element.value_dimension(0)

    # Insert code for computing quantities needed for Piola mapping
    if need_jacobian:
        code.insert(0, format["snippet jacobian"](element.geometric_dimension()) % {"restriction": ""})

    return code
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# From dofmap.py
def __extract_sub_dof_maps(dof_map, parent):
    """Recursively extract sub dof maps as a list of tuples where
    each tuple consists of a tuple labeling the sub dof map and
    the sub dof map itself."""

    if dof_map.num_sub_dof_maps() == 1:
        return [(parent, dof_map)]
    sub_dof_maps = []
    for i in range(dof_map.num_sub_dof_maps()):
        sub_dof_maps += __extract_sub_dof_maps(dof_map.sub_dof_map(i), parent + (i,))
    return sub_dof_maps + [(parent, dof_map)]

def __generate_needs_mesh_entities(dof_map, format):
    "Generate code for needs_mesh_entities"

    # Get total number of dofs per dimension
    num_dofs_per_dim = dof_map.num_dofs_per_dim()

    # Entities needed if at least one dof is associated
    code = [format["bool"](num_dofs_per_dim[dim] > 0) for dim in range(len(num_dofs_per_dim))]

    return code

def __generate_global_dimension(dof_map, format):
    "Generate code for global dimension"

    # Get total number of dofs per dimension
    num_dofs_per_dim = dof_map.num_dofs_per_dim()

    # Sum the number of dofs for each dimension
    terms = []
    for dim in range(len(num_dofs_per_dim)):
        n = num_dofs_per_dim[dim]
        if n == 1:
            terms += [format["num entities"](dim)]
        elif n > 1:
            terms += [format["multiply"]([str(n), format["num entities"](dim)])]

    # Special case, no terms
    if len(terms) == 0:
        code = "0"
    else:
        code = format["add"](terms)

    return code

def __generate_tabulate_dofs(dof_map, format, skip=[]):
    "Generate code for tabulate_dofs"

    # Generate code as a list of declarations
    code = []

    # Iterate over sub dofs
    offset_declared = False
    offset_code = []
    local_offset = 0
    for sub_dof_map in range(len(dof_map.entity_dofs())):

        # Get entity dofs for sub dof map
        sub_entity_dofs = dof_map.entity_dofs()[sub_dof_map]

        # Get the number of dofs per dimension for sub dof map
        num_dofs_per_dim = dof_map.num_dofs_per_dim(sub_dof_map)

        # Iterate over dimensions
        num_dofs = 0
        for dim in sub_entity_dofs:

            # Skip dimension if there are no dofs
            if num_dofs_per_dim[dim] == 0:
                continue

            # Write offset code
            code += offset_code

            # Iterate over entities in dimension
            for entity in sub_entity_dofs[dim]:

                # Iterate over dofs on entity
                for pos in range(len(sub_entity_dofs[dim][entity])):

                    # Get number of current dof
                    dof = sub_entity_dofs[dim][entity][pos]

                    # Assign dof
                    name = format["dofs"](local_offset + dof)
                    if num_dofs_per_dim[dim] > 1:
                        value = format["multiply"](["%d" % num_dofs_per_dim[dim], format["entity index"](dim, entity)])
                    else:
                        value = format["entity index"](dim, entity)

                    # Add position on entity if any
                    if pos > 0:
                        value = format["add"]([value, "%d" % pos])

                    # Add offset if any
                    if offset_declared:
                        value = format["add"]([format["offset access"], value])

                    # Add declaration
                    if (dim, entity) not in skip:
                        code += [(name, value)]

                    # Count the number of dofs for sub dof map
                    num_dofs += 1

            # Update offset
            if num_dofs_per_dim[dim] > 0:

                # Compute additional offset
                if num_dofs_per_dim[dim] > 1:
                    value = format["multiply"](["%d" % num_dofs_per_dim[dim], format["num entities"](dim)])
                else:
                    value = format["num entities"](dim)

                # Add to previous offset
                if not offset_declared:
                    name = format["offset declaration"]
                    offset_declared = True
                else:
                    name = format["offset access"]
                    value = format["add"]([name, value])

                offset_code = [(name, value)]

        # Add to local offset
        local_offset += num_dofs

    return code

def __generate_tabulate_facet_dofs(dof_map, format):
    "Generate code for tabulate_dofs"

    # Get the number of facets
    num_facets = dof_map.element().cell().num_facets()

    # Get incidence
    incidence = dof_map.incidence()

    # Get topological dimension
    D = max([pair[0][0] for pair in incidence])

    # Find out which entities are incident to each facet
    incident = num_facets*[[]]
    for facet in range(num_facets):
        incident[facet] = [pair[1] for pair in incidence if incidence[pair] == True and pair[0] == (D - 1, facet)]

    # Tabulate dofs for each facet
    code = []
    for facet in range(num_facets):
        case = []
        facet_dof = 0
        local_offset = 0
        for sub_entity_dofs in dof_map.entity_dofs():
            num_dofs = 0
            for dim in sub_entity_dofs:
                for entity in sub_entity_dofs[dim]:
                    for dof in sub_entity_dofs[dim][entity]:
                        if (dim, entity) in incident[facet]:
                            name = format["dofs"](facet_dof)
                            value = "%d" % (local_offset + dof)
                            case += [(name, value)]
                            facet_dof += 1
                        num_dofs += 1
            local_offset += num_dofs
        code += [case]

    return code

def __generate_tabulate_coordinates(dof_map, format):
    "Generate code to compute coordinates of the dofs"

    code = []

    # Prefetch formats to speed up code generation
    format_coordinates          = format["argument coordinates"]
    format_element_coordinates  = format["element coordinates"]
    format_matrix_access        = format["matrix access"]

    # Get coordinates of the dofs (on FFC reference element)
    coordinates = dof_map.dof_coordinates()

    # Check if we get some points from fem.dofmap.py
    if not None in coordinates:
        #code += [format["comment"]("This function is implemented assuming affine mapping!!")]
        #code += [format["comment"]("Get cell vertices")]
        code += [format["get cell vertices"]]

        # Create linear Lagrange element for the transformation
        ufl_element = UFLFiniteElement("Lagrange", dof_map.element().cell().domain(), 1)
        element = FiniteElement(ufl_element)

        # Tabulate values of basisfunctions
        table = element.tabulate(0, coordinates)

        # Get the cell shape TODO: KBO: should it be topological_dimension?
        cell_shape = dof_map.element().geometric_dimension()

        # Get matrix of values of basisfunctions at points (dof, values at dofs on linear element)
        transformed_values = numpy.transpose(table[0][(0,)*cell_shape])

        # Get shape of matrix
        shape_val = numpy.shape(transformed_values)

        # Loop dofs
        for i in range(shape_val[0]):
            for j in range(cell_shape):

                name = format_coordinates + format_matrix_access(i,j)

                values = [transformed_values[i][k] for k in range(shape_val[1])]
                symbols = [format_element_coordinates(k,j) for k in range(shape_val[1])]

                # Let inner_product handle the format
                value = inner_product(values, symbols, format)

                code += [(name, value)]

    else:
        code += [format["exception"]("tabulate_coordinates not implemented for this type of element")]

    return code
#------------------------------------------------------------------------------

