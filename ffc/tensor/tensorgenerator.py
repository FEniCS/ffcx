"Code generator for tensor representation"

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2004-11-03"
__copyright__ = "Copyright (C) 2004-2010 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Kristian B. Oelgaard, 2009
# Modified by Marie Rognes (meg@math.uio.no), 2007
# Modified by Garth N. Wells, 2009
# Last changed: 2010-01-13

# FFC modules
from ffc.log import info
from ffc.cpp import format, remove_unused, count_ops
from ffc import codesnippets

# FFC tensor representation modules
from ffc.tensor.monomialtransformation import MonomialIndex

def generate_integral_code(integral_ir,
                           integral_type,
                           sub_domain,
                           ir, prefix, options):
    "Generate code for integral from intermediate representation."

    # Prefetch formatting to speedup code generation
    do_nothing = format["do nothing"]

    # Generate code
    code = {}
    code["classname"] = "%s_%s_%d" % (prefix.lower(), integral_type, sub_domain)
    code["members"] = ""
    code["constructor"] = do_nothing
    code["destructor"] = do_nothing
    code["tabulate_tensor"] = _tabulate_tensor(integral_ir, integral_type, ir, options)

    return code

def _tabulate_tensor(integral_ir, integral_type, ir, options):
    "Generate code for tabulate_tensor."

    # Prefetch formats to speed up code generation
    comment = format["comment"]
    switch = format["switch"]

    # Set of used variables for Jacobian and geometry tensor
    j_set = set()
    g_set = set()

    # Generate code for tensor contraction and geometry tensor
    if integral_type == "cell_integral":

        # Generate code for one single tensor contraction
        t_code = _generate_tensor_contraction(integral_ir, options, g_set)

        # Generate code for geometry tensors
        g_code = _generate_geometry_tensors(integral_ir, j_set, g_set)

    elif integral_type == "interior_facet_integral":

        # Generate code for num_facets tensor contractions
        cases = [None for i in range(ir.num_facets)]
        for i in range(ir.num_facets):
            cases[i] = _generate_tensor_contraction(integral_ir[i], options, g_set)
        t_code = switch("i", cases)

        # Generate code for geometry tensors
        g_code = _generate_geometry_tensors(integral_ir[0], j_set, g_set)

    elif integral_type == "exterior_facet_integral":

        # Generate code for num_facets x num_facets tensor contractions
        cases = [[None for j in range(ir.num_facets)] for i in range(ir.num_facets)]
        for i in range(ir.num_facets):
            for j in range(ir.num_facets):
                cases[i][j] = _generate_tensor_contraction(integral_ir[i][j], options, g_set)
        t_code = switch("i", [switch("j", cases[i]) for i in range(len(cases))])

        # Generate code for geometry tensors
        g_code = _generate_geometry_tensors(integral_ir[0][0], j_set, g_set)

    # Generate code for element tensor, geometry tensor and Jacobian
    j_code = codesnippets.jacobian[ir.geometric_dimension] % {"restriction": ""}

    # Remove unused declarations from Jacobian code
    #jacobi_code = remove_unused(j_code, j_set)

    # FIXME: Missing stuff from old generate_jacobian

    # Compute total number of operations
    j_ops, g_ops, t_ops = [count_ops(c) for c in (j_code, g_code, t_code)]
    total_ops = j_ops + g_ops + t_ops

    # Add generated code
    code = ""
    code += comment("Number of operations (multiply-add pairs) for Jacobian data:      %d\n" % j_ops)
    code += comment("Number of operations (multiply-add pairs) for geometry tensor:    %d\n" % g_ops)
    code += comment("Number of operations (multiply-add pairs) for tensor contraction: %d\n" % t_ops)
    code += comment("Total number of operations (multiply-add pairs):                  %d\n" % total_ops)
    code += j_code
    code += "\n"
    code += comment("Compute geometry tensor\n")
    code += g_code
    code += "\n"
    code += comment("Compute element tensor\n")
    code += t_code

    return code

def _generate_tensor_contraction(terms, options, g_set):
    "Generate code for computation of tensor contraction."

    # Prefetch formats to speed up code generation
    format_add             = format["add"]
    format_iadd            = format["iadd"]
    format_subtract        = format["subtract"]
    format_multiply        = format["multiply"]
    format_element_tensor  = format["element tensor"]
    format_geometry_tensor = format["geometry tensor"]
    format_float           = format["float"]

    # Generate incremental code for now, might be an option later
    incremental = True

    # Get machine precision
    epsilon = options["epsilon"]

    # Get list of primary indices (should be the same so pick first)
    A0, GK = terms[0]
    primary_indices = A0.primary_multi_index.indices

    # Generate code for geometry tensor entries
    gk_tensor = []
    for (j, (A0, GK)) in enumerate(terms):
        gk_tensor_j = []
        for a in A0.secondary_multi_index.indices:
            gk_tensor_j.append((format_geometry_tensor(j, a), a))
        gk_tensor.append((gk_tensor_j, j))

    # Reset variables
    k = 0
    num_dropped = 0
    zero = format_float(0.0)
    code = ""

    # Generate code for computing the element tensor
    for i in primary_indices:
        name = format_element_tensor(k)
        value = None
        for (gka, j) in gk_tensor:
            (A0, GK) = terms[j]
            for (gk, a) in gka:
                a0 = A0.A0[tuple(i + a)]
                if abs(a0) > epsilon:
                    if value and a0 < 0.0:
                        value = format_subtract([value, format_multiply([format_float(-a0), gk])])
                        g_set.add(gk)
                    elif value:
                        value = format_add([value, format_multiply([format_float(a0), gk])])
                        g_set.add(gk)
                    else:
                        value = format_multiply([format_float(a0), gk])
                        g_set.add(gk)
                else:
                    num_dropped += 1
        value = value or zero
        if len(code) > 0:
            code += "\n"
        if incremental:
            code += format_iadd(name, value)
        else:
            code += format_assign(name, value)
        k += 1

    return code

def _generate_geometry_tensors(terms, j_set, g_set):
    "Generate code for computation of geometry tensors."

    # Prefetch formats to speed up code generation
    format_add             = format["add"]
    format_geometry_tensor = format["geometry tensor"]
    format_scale_factor    = format["scale factor"]
    format_declaration     = format["const float declaration"]

    # Iterate over all terms
    code = ""
    offset = 0
    for (i, term) in enumerate(terms):

        # Get secondary indices
        A0, GK = term
        secondary_indices = GK.secondary_multi_index.indices

        # Hack to keep old code generation based on factorization of GK
        # in case we want to reimplement factorization
        GKs = [GK]

        # Iterate over secondary indices
        for a in secondary_indices:

            # Skip code generation if term is not used
            if not format["geometry tensor"](i, a) in g_set:
                continue

            # Compute factorized values
            values = [_generate_entry(GK, a, offset + j, j_set) for (j, GK) in enumerate(GKs)]

            # Sum factorized values
            name = format_geometry_tensor(i, a)
            value = format_add(values)

            # Multiply with determinant factor
            det = GK.determinant
            value = _multiply_value_by_det(value, GK.determinant, len(values) > 1)

            # Add determinant to transformation set
            #!if dets:
            #!    d0 = [format["power"](format["determinant"](det.restriction),
            #!                          det.power) for det in dets]
            #!    jacobi_set.add(format["multiply"](d0))

            # Add code
            code += format_declaration(name, value) + "\n"

        # Add to offset
        offset += len(GKs)

    # Add scale factor
    j_set.add(format_scale_factor)

    return code

def _generate_entry(GK, a, i, j_set):
    "Generate code for the value of a GK entry."

    # Prefetch formats to speed up code generation
    grouping    = format["grouping"]
    add         = format["add"]
    multiply    = format["multiply"]

    # Compute product of factors outside sum
    factors = _extract_factors(GK, a, None, j_set, MonomialIndex.INTERNAL)

    # Compute sum of products of factors inside sum
    terms = [multiply(_extract_factors(GK, a, b, j_set, MonomialIndex.EXTERNAL))
             for b in GK.external_multi_index.indices]

    # Compute product
    if factors:
        entry = multiply(factors + [grouping(add(terms))])
    else:
        entry = add(terms)

    return entry

def _multiply_value_by_det(value, det, is_sum):
    "Generate code for multiplication of value by determinant."
    if not det.power == 0:
        d = [format["power"](format["determinant"](det.restriction), det.power)]
    else:
        d = []
    if value == "1.0":
        v = []
    elif is_sum:
        v = [format["grouping"](value)]
    else:
        v = [value]
    return format["multiply"](d + [format["scale factor"]] + v)

def _extract_factors(GK, a, b, j_set, index_type):
    "Extract factors of given index type in GK entry."

    # Prefetch formats to speed up code generation
    format_coefficient = format["coefficient"]
    format_transform   = format["transform"]

    # List of factors
    factors = []

    # Compute product of coefficients
    for c in GK.coefficients:
        if index_type == c.index.index_type:
            factors.append(format_coefficient(c.number, c.index(secondary=a)))

    # Compute product of transforms
    for t in GK.transforms:
        if index_type in (t.index0.index_type, t.index1.index_type):
            factors.append(format_transform(t.transform_type,
                                            t.index0(secondary=a, external=b),
                                            t.index1(secondary=a, external=b),
                                            t.restriction))
            j_set.add(factors[-1])

    return factors
