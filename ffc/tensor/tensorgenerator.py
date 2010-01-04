"Code generator for tensor representation"

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2004-11-03"
__copyright__ = "Copyright (C) 2004-2009 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Kristian B. Oelgaard, 2009
# Modified by Marie Rognes (meg@math.uio.no), 2007
# Modified by Garth N. Wells, 2009
# Last changed: 2010-01-04

# FFC modules
from ffc.log import info
from ffc.cpp import format, remove_unused
from ffc import codesnippets

# FFC tensor representation modules
from ffc.tensor.monomialtransformation import MonomialIndex

def generate_tensor_integrals(ir, options):
    "Generate code for integrals from intermediate representation."

    # Check if code needs to be generated
    if ir.num_integrals == 0: return {}, {}, {}

    info("Generating code using tensor representation")

    # Generate incremental code for now, might be an option later
    incremental = True

    # Generate code for cell integrals
    code_cell_integrals = []
    for (sub_domain, terms) in enumerate(ir.cell_integrals):
        if len(terms) > 0:
            I = _generate_cell_integral(terms, ir, incremental, format)
            code_cell_integrals.append(I)

    # Generate code for exterior facet integrals
    code_exterior_facet_integrals = []
    for (sub_domain, terms) in enumerate(ir.exterior_facet_integrals):
        if len(terms) > 0:
            I = _generate_exterior_facet_integral(terms, ir, incremental)
            code_exterior_facet_integrals.append(I)

    # Generate code for interior facet integrals
    code_interior_facet_integrals = []
    for (sub_domain, terms) in enumerate(ir.interior_facet_integrals):
        if len(terms) > 0:
            I = _generate_interior_facet_integral(terms, ir, incremental)
            code_interior_facet_integrals.append(I)

    return code_cell_integrals, code_exterior_facet_integrals, code_interior_facet_integrals

def _generate_cell_integral(terms, ir, incremental):
    "Generate code for cell integral."

    # Prefetch formats to speed up code generation
    format_comment = format["comment"]

    # Generate code for everything except tabulate_tensor
    code = {}
    code["classname"] = "FooCellIntegral"
    code["members"] = ""
    code["constructor"] = ""
    code["destructor"] = ""

    code["tabulate_tensor"] = ""

    return code

    # Generate code for element tensor, geometry tensor and Jacobian
    tensor_code, tensor_ops, geometry_set   = _generate_element_tensor(terms, incremental)
    geometry_code, geometry_ops, jacobi_set = _generate_geometry_tensors(terms, geometry_set)
    jacobi_code = codesnippets.jacobian[ir.geometric_dimension]

    # Remove unused declarations from Jacobian code
    jacobi_code = _remove_unused(jacobi_code, jacobi_set)

    # FIXME: Missing stuff from old generate_jacobian

    # Compute total number of operations
    total_ops = tensor_ops + geometry_ops
    code += [format_comment("Number of operations to compute geometry tensor:     %d" % geometry_ops)]
    code += [format_comment("Number of operations to compute tensor contraction:  %d" % tensor_ops)]
    code += [format_comment("Total number of operations to compute cell tensor:   %d" % total_ops)]
    code += [""]
    info("Number of operations to compute tensor: %d" % total_ops)

    # Add generated code
    code += jacobi_code
    code += [format_comment("Compute geometry tensor")]
    code += geometry_code
    code += [""]
    code += [format_comment("Compute element tensor")]
    code += tensor_code

    return {"tabulate_tensor": code, "members": ""}

def _generate_exterior_facet_integral(terms, ir, incremental):
    "Generate code for exterior facet integral."

    # Generate code
    code = {}
    code["classname"] = "FooExteriorFacetIntegral"
    code["members"] = ""
    code["constructor"] = ""
    code["destructor"] = ""
    code["tabulate_tensor"] = ""

    return code

    # Special case: zero contribution
    if all([len(t) == 0 for t in terms]):
        return {"tabulate_tensor": ([format_comment("Do nothing")], []), "members": ""}

    # Generate tensor code + set of used geometry terms
    num_facets = form_representation.num_facets
    cases = [None for i in range(num_facets)]
    geometry_set = set()
    tensor_ops = 0
    for i in range(form_representation.num_facets):
        cases[i], t_ops, g_set = _generate_element_tensor(terms[i], incremental, format)
        geometry_set = geometry_set.union(g_set)
        tensor_ops += t_ops
    tensor_ops = float(tensor_ops) / float(form_representation.num_facets)

    # Generate geometry code + set of used jacobi terms (should be the same, so pick first)
    geometry_code, geometry_ops, jacobi_set = _generate_geometry_tensors(terms[0], geometry_set, format)

    # Generate code for Jacobian
    jacobi_code = [format["generate jacobian"](form_representation.geometric_dimension, "exterior facet")]
    jacobi_code = _remove_unused(jacobi_code, jacobi_set, format)

    # Compute total number of operations
    total_ops = tensor_ops + geometry_ops
    code += [format_comment("Number of operations to compute geometry tensor:     %d" % geometry_ops)]
    code += [format_comment("Number of operations to compute tensor contraction:  %d" % tensor_ops)]
    code += [format_comment("Total number of operations to compute facet tensor:  %d" % total_ops)]
    code += [""]
    info("Number of operations to compute tensor: %d" % total_ops)

    # Add generated code
    code += jacobi_code
    code += [format_comment("Compute geometry tensor")]
    code += geometry_code
    code += [""]
    code += [format_comment("Compute element tensor")]

    return {"tabulate_tensor": (code, cases), "members": ""}

def _generate_interior_facet_integral(terms, ir, incremental):
    "Generate code for interior facet integral."

    # Generate code
    code = {}
    code["classname"] = "FooInteriorFacetIntegral"
    code["members"] = ""
    code["constructor"] = ""
    code["destructor"] = ""
    code["tabulate_tensor"] = ""

    return code

    # Special case: zero contribution
    if all([len(t) == 0 for tt in terms for t in tt]):
        return {"tabulate_tensor": ([format_comment("Do nothing")], []), "members": ""}

    # Special case: zero contribution
    if len(terms) == 0: return {"tabulate_tensor": "", "members": ""}

    # Generate tensor code + set of used geometry terms
    num_facets = form_representation.num_facets
    cases = [[None for j in range(num_facets)] for i in range(num_facets)]
    geometry_set = set()
    tensor_ops = 0
    for i in range(num_facets):
        for j in range(num_facets):
            cases[i][j], t_ops, g_set = _generate_element_tensor(terms[i][j], incremental, format)
            geometry_set = geometry_set.union(g_set)
            tensor_ops += t_ops
    tensor_ops = float(tensor_ops) / float(form_representation.num_facets)

    # Generate geometry code + set of used jacobi terms (should be the same, so pick first)
    geometry_code, geometry_ops, jacobi_set = _generate_geometry_tensors(terms[0][0], geometry_set, format)

    # Generate code for Jacobian
    jacobi_code = [format["generate jacobian"](form_representation.geometric_dimension, "interior facet")]
    jacobi_code = _remove_unused(jacobi_code, jacobi_set, format)

    # Compute total number of operations
    total_ops = tensor_ops + geometry_ops
    code += [format_comment("Number of operations to compute geometry tensor:     %d" % geometry_ops)]
    code += [format_comment("Number of operations to compute tensor contraction:  %d" % tensor_ops)]
    code += [format_comment("Total number of operations to compute facet tensor:  %d" % total_ops)]
    code += [""]
    info("Number of operations to compute tensor: %d" % total_ops)

    # Add generated code
    code += jacobi_code
    code += [format_comment("Compute geometry tensor")]
    code += geometry_code
    code += [""]
    code += [format_comment("Compute element tensor")]

    return {"tabulate_tensor": (code, cases), "members": ""}

def _generate_element_tensor(terms, incremental, format):
    "Generate list of declarations for computation of element tensor."

    # Generate code as a list of declarations
    code = []

    # Prefetch formats to speed up code generation
    format_element_tensor  = format["element tensor"]
    format_geometry_tensor = format["geometry tensor"]
    format_add             = format["add"]
    format_iadd            = format["iadd"]
    format_subtract        = format["subtract"]
    format_multiply        = format["multiply"]
    format_float           = format["float"]
    format_epsilon         = format["epsilon"]

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

    # Generate code for computing the element tensor
    k = 0
    num_dropped = 0
    num_ops = 0
    zero = format_float(0.0)
    geometry_set = set()
    for i in primary_indices:
        name = format_element_tensor(k)
        value = None
        for (gka, j) in gk_tensor:
            (A0, GK) = terms[j]
            for (gk, a) in gka:
                a0 = A0.A0[tuple(i + a)]
                if abs(a0) > format_epsilon:
                    if value and a0 < 0.0:
                        value = format_subtract([value, format_multiply([format_floating_point(-a0), gk])])
                        geometry_set.add(gk)
                        num_ops += 1
                    elif value:
                        value = format_add([value, format_multiply([format_floating_point(a0), gk])])
                        geometry_set.add(gk)
                        num_ops += 1
                    else:
                        value = format_multiply([format_floating_point(a0), gk])
                        geometry_set.add(gk)
                    num_ops += 1
                else:
                    num_dropped += 1
        value = value or zero
        if incremental:
            code += [format_iadd(name, value)]
        else:
            code += [(name, value)]
        k += 1

    return (code, num_ops, geometry_set)

def _generate_geometry_tensors(terms, geometry_set, format):
    "Generate list of declarations for computation of geometry tensors"

    # Generate code as a list of declarations
    code = []

    # Prefetch formats to speed up code generation
    format_add             = format["add"]
    format_geometry_tensor = format["geometry tensor"]
    format_scale_factor    = format["scale factor"]

    # Iterate over all terms
    j = 0
    jacobi_set = set()
    num_ops = 0
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
                val, entry_ops, t_set = _generate_entry(GK, a, jj, format)
                values += [val]
                num_ops += entry_ops
                jacobi_set = jacobi_set | t_set
                jj += 1

            # Sum factorized values
            if values:
                num_ops += len(values) - 1
            name = format_geometry_tensor(i, a) # FIXME: Need declaration here
            value = format_add(values)

            # Multiply with determinant factor
            det = GK.determinant
            value = _multiply_value_by_det(value, GK.determinant, format, len(values) > 1)
            num_ops += 1

            # Add determinant to transformation set
            #!if dets:
            #!    d0 = [format["power"](format["determinant"](det.restriction),
            #!                          det.power) for det in dets]
            #!    jacobi_set.add(format["multiply"](d0))

            # Add declaration
            code += [(name, value)]

        j += len(GKs)

    # Add scale factor
    jacobi_set.add(format_scale_factor)

    return (code, num_ops, jacobi_set)

def _generate_entry(GK, a, i, format):
    "Generate code for the value of entry a of geometry tensor G"

    jacobi_set = set()

    # Compute product of factors outside sum
    factors = []
    num_ops = 0
    for j in range(len(GK.coefficients)):
        c = GK.coefficients[j]
        if not c.index.index_type == MonomialIndex.EXTERNAL:
            #coefficient = format["coefficient"](c.number, i, j, c.index(secondary=a))
            coefficient = format["coefficient"](c.number, c.index(secondary=a))
            factors += [coefficient]

    for t in GK.transforms:
        if not (t.index0.index_type == MonomialIndex.EXTERNAL or t.index1.index_type == MonomialIndex.EXTERNAL):
            trans = format["transform_ufl"](t.transform_type,
                                            t.index0(secondary=a),
                                            t.index1(secondary=a),
                                            t.restriction)
            factors += [trans]
            jacobi_set.add(trans)

    if factors:
        num_ops += len(factors) - 1

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
                coefficient = format["coefficient"](c.number, c.index(secondary=a))
                factors += [coefficient]
        for t in GK.transforms:
            if t.index0.index_type == MonomialIndex.EXTERNAL or t.index1.index_type == MonomialIndex.EXTERNAL:
                trans = format["transform_ufl"](t.transform_type,
                                                t.index0(secondary=a, external=b),
                                                t.index1(secondary=a, external=b),
                                                t.restriction)
                factors += [trans]
                jacobi_set.add(trans)
        if factors:
            num_ops += len(factors) - 1
        terms += [format["multiply"](factors)]

    if terms:
        num_ops += len(terms) - 1

    sum = format["add"](terms)
    if sum: sum = format["grouping"](sum)
    if sum: f1 = [sum]
    else: f1 = []

    fs = f0 + f1
    if not fs: fs = ["1.0"]
    else:
        num_ops += len(fs) - 1

    # Compute product of all factors
    return (format["multiply"](fs), num_ops, jacobi_set)

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

def _remove_unused(code, set):
    """
    Remove unused variables so that the compiler will not complain.

    Normally, the removal of unused variables should happen at the
    formatting stage, but since the code for the tensor contraction
    may grow to considerable size, we make an exception and remove
    unused variables here when we know the names of the used
    variables. No searching necessary and much, much, much faster.
    """

    return code

    if not code: return code

    # Generate body of code, using the format
    lines = format["generate body"](code)

    # Generate auxiliary code line that uses all members of the set (to trick remove_unused)
    line_set = format["iadd"]("A", format["multiply"](set))
    lines += "\n" + line_set

    # Remove unused Jacobi declarations
    code = remove_unused(lines)

    # Delete auxiliary line
    code = code.replace("\n" + line_set, "")

    return [code]
