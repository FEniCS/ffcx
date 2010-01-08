"Code generator for tensor representation"

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2004-11-03"
__copyright__ = "Copyright (C) 2004-2010 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Kristian B. Oelgaard, 2009
# Modified by Marie Rognes (meg@math.uio.no), 2007
# Modified by Garth N. Wells, 2009
# Last changed: 2010-01-08

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

    g_set = set()

    # Generate code for tensor contraction and geometry tensor
    if integral_type == "cell_integral":

        # Generate code for one single tensor contraction
        t_code = _generate_tensor_contraction(integral_ir, options, g_set)

        # Generate code for geometry tensors
        g_code, j_set = _generate_geometry_tensors(integral_ir, g_set)

    elif integral_type == "interior_facet_integral":

        # Generate code for num_facets tensor contractions
        cases = [None for i in range(ir.num_facets)]
        for i in range(ir.num_facets):
            cases[i] = _generate_tensor_contraction(integral_ir[i], options, g_set)
        t_code = switch("i", cases)

        # Generate code for geometry tensors
        g_code, j_set = _generate_geometry_tensors(integral_ir[0], g_set)

    elif integral_type == "exterior_facet_integral":

        # Generate code for num_facets x num_facets tensor contractions
        cases = [[None for j in range(ir.num_facets)] for i in range(ir.num_facets)]
        for i in range(ir.num_facets):
            for j in range(ir.num_facets):
                cases[i][j] = _generate_tensor_contraction(integral_ir[i][j], options, g_set)
        t_code = switch("i", [switch("j", cases[i]) for i in range(len(cases))])

        # Generate code for geometry tensors
        g_code, j_set = _generate_geometry_tensors(integral_ir[0][0], g_set)

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

def _generate_geometry_tensors(terms, geometry_set):
    "Generate code for computation of geometry tensors."

    # Prefetch formats to speed up code generation
    format_add             = format["add"]
    format_geometry_tensor = format["geometry tensor"]
    format_scale_factor    = format["scale factor"]
    format_declaration     = format["const float declaration"]

    # Reset variables
    j = 0
    jacobi_set = set()
    code = ""

    # Iterate over all terms
    for (i, term) in enumerate(terms):

        # Get list of secondary indices (should be the same so pick first)
        A0, GK = term
        secondary_indices = GK.secondary_multi_index.indices

        # Hack to keep old code generation based on factorization of GK
        # in case we want to reimplement factorization
        GKs = [GK]

        # Iterate over secondary indices
        for a in secondary_indices:

            # Skip code generation if term is not used
            if not format["geometry tensor"](i, a) in geometry_set:
                continue

            # Compute factorized values
            values = []
            jj = j
            for GK in GKs:
                val, t_set = _generate_entry(GK, a, jj, format)
                values += [val]
                jacobi_set = jacobi_set | t_set
                jj += 1

            # Sum factorized values
            name = format_geometry_tensor(i, a) # FIXME: Need declaration here
            value = format_add(values)

            # Multiply with determinant factor
            det = GK.determinant
            value = _multiply_value_by_det(value, GK.determinant, format, len(values) > 1)

            # Add determinant to transformation set
            #!if dets:
            #!    d0 = [format["power"](format["determinant"](det.restriction),
            #!                          det.power) for det in dets]
            #!    jacobi_set.add(format["multiply"](d0))

            # Add code
            code += format_declaration(name, value) + "\n"

        j += len(GKs)

    # Add scale factor
    jacobi_set.add(format_scale_factor)

    return (code, jacobi_set)

def _generate_entry(GK, a, i, format):
    "Generate code for the value of entry a of geometry tensor G."

    # Prefetch formats to speed up code generation
    format_grouping    = format["grouping"]
    format_add         = format["add"]
    format_multiply    = format["multiply"]
    format_coefficient = format["coefficient"]
    format_transform   = format["transform"]

    # Reset variables
    factors = []
    jacobi_set = set()

    # Compute product of factors outside sum
    for j in range(len(GK.coefficients)):
        c = GK.coefficients[j]
        if not c.index.index_type == MonomialIndex.EXTERNAL:
            coefficient = format_coefficient(c.number, c.index(secondary=a))
            factors += [coefficient]

    for t in GK.transforms:
        if not (t.index0.index_type == MonomialIndex.EXTERNAL or t.index1.index_type == MonomialIndex.EXTERNAL):
            trans = format_transform(t.transform_type,
                                     t.index0(secondary=a),
                                     t.index1(secondary=a),
                                     t.restriction)
            factors += [trans]
            jacobi_set.add(trans)

    monomial = format["multiply"](factors)
    if monomial: f0 = [monomial]
    else: f0 = []

    # Compute sum of monomials inside sum
    terms = []
    for b in GK.external_multi_index.indices:
        factors = []
        for j in range(len(GK.coefficients)):
            c = GK.coefficients[j]
            if c.index.index_type == MonomialIndex.EXTERNAL:
                coefficient = format_coefficient(c.number, c.index(secondary=a))
                factors += [coefficient]
        for t in GK.transforms:
            if t.index0.index_type == MonomialIndex.EXTERNAL or t.index1.index_type == MonomialIndex.EXTERNAL:
                trans = format_transform(t.transform_type,
                                         t.index0(secondary=a, external=b),
                                         t.index1(secondary=a, external=b),
                                         t.restriction)
                factors += [trans]
                jacobi_set.add(trans)
        terms += [format_multiply(factors)]

    sum = format_add(terms)
    if sum: sum = format_grouping(sum)
    if sum: f1 = [sum]
    else: f1 = []

    fs = f0 + f1
    if not fs: fs = ["1.0"]

    # Compute product of all factors
    return (format_multiply(fs), jacobi_set)

def _multiply_value_by_det(value, det, format, is_sum):
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
