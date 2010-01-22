"Code generator for tensor representation"

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2004-11-03"
__copyright__ = "Copyright (C) 2004-2010 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Kristian B. Oelgaard, 2009
# Modified by Marie Rognes (meg@math.uio.no), 2007
# Modified by Garth N. Wells, 2009
# Last changed: 2010-01-19

# FFC modules
from ffc.cpp import format, remove_unused, count_ops
from ffc import codesnippets

# FFC tensor representation modules
from ffc.tensor.monomialtransformation import MonomialIndex

def generate_integral_code(ir, prefix, options):
    "Generate code for integral from intermediate representation."

    # Prefetch formatting to speedup code generation
    do_nothing = format["do nothing"]
    classname = format["classname " + ir["integral_type"]]

    # Generate code
    code = {}
    code["classname"] = classname(prefix, ir["form_id"], ir["sub_domain"])
    code["members"] = ""
    code["constructor"] = do_nothing
    code["destructor"] = do_nothing
    code["tabulate_tensor"] = _tabulate_tensor(ir, options)

    return code

def _tabulate_tensor(ir, options):
    "Generate code for tabulate_tensor."

    # Prefetch formats to speed up code generation
    comment = format["comment"]
    switch = format["switch"]

    # Set of used variables for Jacobian and geometry tensor
    j_set = set()
    g_set = set()

    # Extract data from intermediate representation
    AK = ir["AK"]
    integral_type = ir["integral_type"]
    geometric_dimension = ir["geometric_dimension"]
    num_facets = ir["num_facets"]

    # Check integral typpe and generate code
    if integral_type == "cell_integral":

        # Generate code for one single tensor contraction
        t_code = _generate_tensor_contraction(AK, options, g_set)

        # Generate code for geometry tensors
        g_code = _generate_geometry_tensors(AK, j_set, g_set)

        # Generate code for Jacobian
        r = {"restriction": ""}
        j_code  = codesnippets.jacobian[geometric_dimension] % r
        j_code += "\n\n" + codesnippets.scale_factor

    elif integral_type == "exterior_facet_integral":

        # Generate code for num_facets tensor contractions
        cases = [None for i in range(num_facets)]
        for i in range(num_facets):
            cases[i] = _generate_tensor_contraction(AK[i], options, g_set)
        t_code = switch("facet", cases)

        # Generate code for geometry tensors
        g_code = _generate_geometry_tensors(AK[0], j_set, g_set)

        # Generate code for Jacobian
        r = {"restriction": ""}
        j_code  = codesnippets.jacobian[geometric_dimension] % r
        j_code += "\n\n" + codesnippets.facet_determinant[geometric_dimension] % r

    elif integral_type == "interior_facet_integral":

        # Generate code for num_facets x num_facets tensor contractions
        cases = [[None for j in range(num_facets)] for i in range(num_facets)]
        for i in range(num_facets):
            for j in range(num_facets):
                cases[i][j] = _generate_tensor_contraction(AK[i][j], options, g_set)
        t_code = switch("facet0", [switch("facet1", cases[i]) for i in range(len(cases))])

        # Generate code for geometry tensors
        g_code = _generate_geometry_tensors(AK[0][0], j_set, g_set)

        # Generate code for Jacobian
        r = {"restriction": ""}
        j_code  = codesnippets.jacobian[ir.geometric_dimension] % r
        j_code += "\n\n" + codesnippets.facet_determinant[ir.geometric_dimension] % r

    else:
        error("Unhandled integral type: " + str(integral_type))

    # Remove unused declarations from Jacobian code
    jacobi_code = remove_unused(j_code, j_set)

    # FIXME: Missing stuff from old generate_jacobian

    # Compute total number of operations
    j_ops, g_ops, t_ops = [count_ops(c) for c in (j_code, g_code, t_code)]
    total_ops = j_ops + g_ops + t_ops

    # Add generated code
    lines = []
    lines.append(comment("Number of operations (multiply-add pairs) for Jacobian data:      %d" % j_ops))
    lines.append(comment("Number of operations (multiply-add pairs) for geometry tensor:    %d" % g_ops))
    lines.append(comment("Number of operations (multiply-add pairs) for tensor contraction: %d" % t_ops))
    lines.append(comment("Total number of operations (multiply-add pairs):                  %d" % total_ops))
    lines.append("")
    lines.append(jacobi_code)
    lines.append("")
    lines.append(comment("Compute geometry tensor"))
    lines.append(g_code)
    lines.append("")
    lines.append(comment("Compute element tensor"))
    lines.append(t_code)

    return "\n".join(lines)

def _generate_tensor_contraction(terms, options, g_set):
    "Generate code for computation of tensor contraction."

    # Prefetch formats to speed up code generation
    add             = format["addition"]
    iadd            = format["iadd"]
    subtract        = format["sub"]
    multiply        = format["multiply"]
    assign          = format["assign"]
    element_tensor  = format["element tensor"]
    geometry_tensor = format["geometry tensor"]
    format_float    = format["float"]
    zero            = format_float(0)

    # FIXME: This should be an option
    incremental = False

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
            gk_tensor_j.append((geometry_tensor(j, a), a))
        gk_tensor.append((gk_tensor_j, j))

    # Generate code for computing the element tensor
    lines = []
    for (k, i) in enumerate(primary_indices):
        name = element_tensor(k)
        value = None
        for (gka, j) in gk_tensor:
            (A0, GK) = terms[j]
            for (gk, a) in gka:
                a0 = A0.A0[tuple(i + a)]

                # Skip small values
                if abs(a0) < epsilon: continue

                # Compute value
                if value and a0 < 0.0:
                    value = subtract([value, multiply([format_float(-a0), gk])])
                elif value:
                    value = add([value, multiply([format_float(a0), gk])])
                else:
                    value = multiply([format_float(a0), gk])

                # Remember that gk has been used
                g_set.add(gk)

        # Handle special case
        value = value or zero

        # Add value
        if incremental:
            lines.append(iadd(name, value))
        else:
            lines.append(assign(name, value))

    return "\n".join(lines)

def _generate_geometry_tensors(terms, j_set, g_set):
    "Generate code for computation of geometry tensors."

    # Prefetch formats to speed up code generation
    format_add             = format["addition"]
    format_geometry_tensor = format["geometry tensor"]
    format_scale_factor    = format["scale factor"]
    format_declaration     = format["const float declaration"]

    # Iterate over all terms
    lines = []
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
            if not format["geometry tensor"](i, a) in g_set: continue

            # Compute factorized values
            values = [_generate_entry(GK, a, offset + j, j_set) for (j, GK) in enumerate(GKs)]

            # Sum factorized values
            name = format_geometry_tensor(i, a)
            value = format_add(values)

            # Multiply with determinant factor
            det = GK.determinant
            value = _multiply_value_by_det(value, GK.determinant, len(values) > 1)

            # Add code
            lines.append(format_declaration(name, value))

        # Add to offset
        offset += len(GKs)

    # Add scale factor
    j_set.add(format_scale_factor)

    return "\n".join(lines)

def _generate_entry(GK, a, i, j_set):
    "Generate code for the value of a GK entry."

    # Prefetch formats to speed up code generation
    grouping = format["grouping"]
    add      = format["addition"]
    multiply = format["multiply"]

    # Compute product of factors outside sum
    factors = _extract_factors(GK, a, None, j_set, index_type=MonomialIndex.SECONDARY)

    # Compute sum of products of factors inside sum
    terms = [multiply(_extract_factors(GK, a, b, j_set, index_type=MonomialIndex.EXTERNAL))
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
        d = [format["power"](format["det(J)"](det.restriction), det.power)]
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
    coefficient = format["coefficient"]
    transform   = format["transform"]

    # List of factors
    factors = []

    # Compute product of coefficients
    for c in GK.coefficients:
        if c.index.index_type == index_type:
            factors.append(coefficient(c.number, c.index(secondary=a)))

    # Compute product of transforms
    for t in GK.transforms:

        # Note non-trivial logic here
        if index_type == MonomialIndex.EXTERNAL:
            include_index = MonomialIndex.EXTERNAL in (t.index0.index_type, t.index1.index_type)
        else:
            include_index = not (MonomialIndex.EXTERNAL in (t.index0.index_type, t.index1.index_type))

        # Add factor
        if include_index:
            factors.append(transform(t.transform_type,
                                     t.index0(secondary=a, external=b),
                                     t.index1(secondary=a, external=b),
                                     t.restriction))
            j_set.add(factors[-1])

    return factors
