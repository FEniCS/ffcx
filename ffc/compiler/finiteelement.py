"Code generation for finite element"

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-01-23 -- 2009-08-26"
__copyright__ = "Copyright (C) 2007-2009 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Kristian Oelgaard 2009
# Modified by Marie Rognes 2008
# Modified by Garth N. Wells 2009

# FIAT modules
from FIAT.shapes import LINE

# FFC common modules
from ffc.common.log import debug, error

# FFC fem modules
from ffc.fem.finiteelement import *
from ffc.fem.vectorelement import *
from ffc.fem.dofmap import *
from ffc.fem.quadratureelement import *

# FFC code generation common modules
from evaluatebasis import evaluate_basis
from evaluatebasisderivatives import evaluate_basis_derivatives
from codeutils import inner_product

def generate_finite_elements(form_data, format):
    "Generate code for finite elements, including recursively nested sub elements."

    debug("Generating code for finite elements...")
    code = []

    # Iterate over form elements
    for (i, element) in enumerate(form_data.ffc_elements):

        # Extract sub elements
        sub_elements = _extract_sub_elements(element, (i,))

        # Generate code for each element
        for (label, sub_element) in sub_elements:
            code += [(label, generate_finite_element(sub_element, format))]
                
    debug("done")
    return code

def generate_finite_element(element, format):
    """Generate dictionary of code for the given finite element
    according to the given format."""

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
    if not True in [isinstance(e, QuadratureElement) for e in element.extract_elements()]:
        # Generate code for evaluate_basis
        code["evaluate_basis"] = evaluate_basis(element, format)

        # Generate code for evaluate_basis_derivatives
        code["evaluate_basis_derivatives"] = evaluate_basis_derivatives(element, format)

        # Generate vectorised version of evaluate functions
        code["evaluate_basis_all"] =\
          format["exception"]("The vectorised version of evaluate_basis() is not yet implemented.")
        code["evaluate_basis_derivatives_all"] =\
          format["exception"]("The vectorised version of evaluate_basis_derivatives() is not yet implemented.")

        # Generate code for interpolate_vertex_values
        code["interpolate_vertex_values"] = __generate_interpolate_vertex_values(element, format)
    else:
        code["evaluate_basis"] =\
          format["exception"]("evaluate_basis() is not supported for QuadratureElement")
        code["evaluate_basis_derivatives"] =\
          format["exception"]("evaluate_basis_derivatives() is not supported for QuadratureElement")

        # Generate vectorised version of evaluate functions
        code["evaluate_basis_all"] =\
          format["exception"]("evaluate_basis_all() is not supported for QuadratureElement.")
        code["evaluate_basis_derivatives_all"] =\
          format["exception"]("evaluate_basis_derivatives_all() is not supported for QuadratureElement.")

        code["interpolate_vertex_values"] =\
          format["exception"]("interpolate_vertex_values() is not supported for QuadratureElement")

    # Generate code for evaluate_dof
    code["evaluate_dof"] = __generate_evaluate_dof(element, format)

    # Generate code for num_sub_elements
    code["num_sub_elements"] = "%d" % element.num_sub_elements()

    return code

def _extract_sub_elements(element, parent):
    """Recursively extract sub elements as a list of tuples where
    each tuple consists of a tuple labeling the sub element and
    the sub element itself."""
    
    if element.num_sub_elements() == 1:
        return [(parent, element)]
    sub_elements = []
    for i in range(element.num_sub_elements()):
        sub_elements += _extract_sub_elements(element.sub_element(i), parent + (i,))
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
                                         element.cell_dimension(), s)]
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
    code += [indent(format["snippet map_onto_physical"](element.cell_dimension())
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
    mappings = [mapping_to_int[element.space_mapping(i)]
                for i in range(element.space_dimension())]
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
        precode += [format["snippet jacobian"](element.cell_dimension())
                 % {"restriction":""}]
        precode += ["\ndouble copyofvalues[%d];" % num_values]
    elif copiola_present:
        # If covariant piola: Will need J only
        precode += [format["snippet only jacobian"](element.cell_dimension())
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
    n = element.cell_dimension()
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

    # Set vertices (note that we need to use the FIAT reference cells)
    if element.cell_shape() == LINE:
        vertices = [(0,), (1,)]
    elif element.cell_shape() == TRIANGLE:
        vertices = [(0, 0), (1, 0), (0, 1)]
    elif element.cell_shape() == TETRAHEDRON:
        vertices =  [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]

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
                    coefficients = table[0][sub_element.cell_dimension()*(0,)][:, v]
                    dof_values = [format["dof values"](offset_dof_values + n) for n in range(len(coefficients))]
                    name = format["vertex values"](size*v + offset_vertex_values)
                    value = inner_product(coefficients, dof_values, format)
                    code += [(name, value)]
            else:
                for dim in range(sub_element.value_dimension(0)):
                    for v in range(len(vertices)):
                        coefficients = table[dim][0][sub_element.cell_dimension()*(0,)][:, v]
                        dof_values = [format["dof values"](offset_dof_values + n) for n in range(len(coefficients))]
                        name = format["vertex values"](size*v + offset_vertex_values + dim)
                        value = inner_product(coefficients, dof_values, format)
                        code += [(name, value)]

        elif (mapping == CONTRAVARIANT_PIOLA or mapping == COVARIANT_PIOLA):

            code += [format["comment"]("Evaluate at vertices and use Piola mapping")]

            # Remember to add code later for Jacobian
            need_jacobian = True

            # Check that dimension matches for Piola transform
            if not sub_element.value_dimension(0) == sub_element.cell_dimension():
                error("Vector dimension of basis function does not match for Piola transform.")

            # Get entities for the dofs
            dof_entities = DofMap(sub_element).dof_entities()

            for dim in range(sub_element.value_dimension(0)):
                for v in range(len(vertices)):
                    terms = []
                    for n in range(sub_element.space_dimension()):
                        # Get basis function values at vertices
                        coefficients = [table[j][0][sub_element.cell_dimension()*(0,)][n, v] for j in range(sub_element.value_dimension(0))]

                        if mapping == COVARIANT_PIOLA:
                            # Get row of inverse transpose Jacobian
                            jacobian_row = [format["transform"]("JINV", j, dim, None) for j in range(sub_element.cell_dimension())]
                        else:
                            # mapping == CONTRAVARIANT_PIOLA:
                            # Get row of Jacobian
                            jacobian_row = [format["transform"]("J", j, dim, None) for j in range(sub_element.cell_dimension())]
                            
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
        code.insert(0, format["snippet jacobian"](element.cell_dimension()) % {"restriction": ""})        
    
    return code
