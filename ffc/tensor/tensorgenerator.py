"Code generator for tensor representation"

# Copyright (C) 2004-2013 Anders Logg
#
# This file is part of FFC.
#
# FFC is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# FFC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with FFC. If not, see <http://www.gnu.org/licenses/>.
#
# Modified by Kristian B. Oelgaard, 2009-2010
# Modified by Marie Rognes, 2007
# Modified by Garth N. Wells, 2009
# Modified by Mehdi Nikbakht, 2010
# Modified by Martin Alnaes, 2013
#
# First added:  2004-11-03
# Last changed: 2013-02-10

# FFC modules
from ffc.log import error
from ffc.cpp import format, remove_unused, count_ops

# FFC tensor representation modules
from ffc.tensor.monomialtransformation import MonomialIndex
from ffc.representationutils import initialize_integral_code

def generate_integral_code(ir, prefix, parameters):
    "Generate code for integral from intermediate representation."
    code = initialize_integral_code(ir, prefix, parameters)
    code["tabulate_tensor"] = _tabulate_tensor(ir, parameters)
    return code

def _tabulate_tensor(ir, parameters):
    "Generate code for tabulate_tensor."

    # Prefetch formats to speed up code generation
    comment = format["comment"]
    switch = format["switch"]

    # Set of used variables for Jacobian and geometry tensor
    j_set = set()
    g_set = set()

    # Extract data from intermediate representation
    AK = ir["AK"]
    domain_type = ir["domain_type"]
    tdim = ir["topological_dimension"]
    gdim = ir["geometric_dimension"]
    oriented = ir["needs_oriented"]
    num_facets = ir["num_facets"]

    # Check integral type and generate code
    if domain_type == "cell":

        # Generate code for one single tensor contraction
        t_code = _generate_tensor_contraction(AK, parameters, g_set)

        # Generate code for geometry tensors
        g_code = _generate_geometry_tensors(AK, j_set, g_set, tdim, gdim)

        # Generate code for basic geometric quantities
        j_code  = ""
        j_code += format["compute_jacobian"](tdim, gdim)
        j_code += "\n"
        j_code += format["compute_jacobian_inverse"](tdim, gdim)
        if oriented:
            j_code += format["orientation"](tdim, gdim)
        j_code += "\n"
        j_code += format["scale factor snippet"]

    elif domain_type == "exterior_facet":

        # Generate code for num_facets tensor contractions
        cases = [None for i in range(num_facets)]
        for i in range(num_facets):
            cases[i] = _generate_tensor_contraction(AK[i], parameters, g_set)
        t_code = switch(format["facet"](None), cases)

        # Generate code for geometry tensors
        g_code = _generate_geometry_tensors(AK[0], j_set, g_set, tdim, gdim)

        # Generate code for Jacobian
        j_code = ""
        j_code += format["compute_jacobian"](tdim, gdim)
        j_code += "\n"
        j_code += format["compute_jacobian_inverse"](tdim, gdim)
        if oriented:
            j_code += format["orientation"](tdim, gdim)
        j_code += "\n"
        j_code += format["facet determinant"](tdim, gdim)

    elif domain_type == "interior_facet":

        # Generate code for num_facets x num_facets tensor contractions
        cases = [[None for j in range(num_facets)] for i in range(num_facets)]
        for i in range(num_facets):
            for j in range(num_facets):
                cases[i][j] = _generate_tensor_contraction(AK[i][j], parameters, g_set)
        t_code = switch(format["facet"]("+"), [switch(format["facet"]("-"), cases[i]) for i in range(len(cases))])

        # Generate code for geometry tensors
        g_code = _generate_geometry_tensors(AK[0][0], j_set, g_set, tdim, gdim)

        # Generate code for Jacobian
        j_code = ""
        for _r in ["+", "-"]:
            j_code += format["compute_jacobian"](tdim, gdim, r=_r)
            j_code += "\n"
            j_code += format["compute_jacobian_inverse"](tdim, gdim, r=_r)
            j_code += "\n"
        j_code += format["facet determinant"](tdim, gdim, r="+")
        j_code += "\n"

    else:
        error("Unhandled integral type: " + str(domain_type))

    # Remove unused declarations from Jacobian code
    j_code = remove_unused(j_code, j_set)

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
    lines.append(j_code)
    lines.append("")
    lines.append(comment("Compute geometry tensor"))
    lines.append(g_code)
    lines.append("")
    lines.append(comment("Compute element tensor"))
    lines.append(t_code)

    return "\n".join(lines)

def _generate_tensor_contraction(terms, parameters, g_set):
    """
    Generate code for computation of tensor contraction, choosing
    either standard or optimized contraction.
    """

    # Only check first term, assuming either non or all are optimized
    A0, GK, optimized_contraction = terms[0]
    if optimized_contraction is None:
        return _generate_tensor_contraction_standard(terms, parameters, g_set)
    else:
        return _generate_tensor_contraction_optimized(terms, parameters, g_set)

def _generate_tensor_contraction_standard(terms, parameters, g_set):
    """
    Generate code for computation of tensor contraction using full
    tensor contraction.
    """

    # Prefetch formats to speed up code generation
    iadd            = format["iadd"]
    assign          = format["assign"]
    element_tensor  = format["element tensor"]
    geometry_tensor = format["geometry tensor"]
    zero            = format["float"](0)
    inner_product   = format["inner product"]

    # True if we should add to element tensor (not used)
    incremental = False

    # Get machine precision
    epsilon = parameters["epsilon"]

    # Get list of primary indices (should be the same so pick first)
    A0, GK, optimized_contraction = terms[0]
    primary_indices = A0.primary_multi_index.indices

    # Generate code for geometry tensor entries
    gk_tensor = []
    for (j, (A0, GK, optimized_contraction)) in enumerate(terms):
        gk_tensor_j = []
        for a in A0.secondary_multi_index.indices:
            gk_tensor_j.append((geometry_tensor(j, a), a))
        gk_tensor.append((gk_tensor_j, j))

    # Generate code for computing the element tensor
    lines = []
    for (k, i) in enumerate(primary_indices):
        name = element_tensor(k)
        coefficients = []
        entries = []
        for (gka, j) in gk_tensor:
            (A0, GK, optimized_contraction) = terms[j]
            for (gk, a) in gka:
                a0 = A0.A0[tuple(i + a)]

                # Skip small values
                if abs(a0) < epsilon: continue

                # Compute value
                coefficients.append(a0)
                entries.append(gk)

                # Remember that gk has been used
                g_set.add(gk)

        # Compute inner product
        value = inner_product(coefficients, entries)

        # Handle special case
        value = value or zero

        # Add value
        if incremental:
            lines.append(iadd(name, value))
        else:
            lines.append(assign(name, value))

    return "\n".join(lines)

def _generate_tensor_contraction_optimized(terms, parameters, g_set):
    """
    Generate code for computation of tensor contraction using
    optimized tensor contraction.
    """

    # Prefetch formats to speed up code generation
    assign            = format["assign"]
    geometry_tensor   = format["geometry tensor"]
    inner_product     = format["inner product"]
    float_declaration = format["float declaration"]

    # Handle naming of entries depending on the number of terms
    if len(terms) == 1:
        element_tensor = lambda i, j: format["element tensor"](i)
    else:
        element_tensor = lambda i, j: format["element tensor term"](i, j)

    # Iterate over terms
    lines = []
    for (j, (A0, GK, optimized_contraction)) in enumerate(terms):

        # Check that an optimized contraction has been computed
        if optimized_contraction is None:
            error("Missing optimized tensor contraction.")

        # Declare array if necessary
        if len(terms) > 1:
            num_entries = len(optimized_contraction)
            lines.append("%s %s;" % (float_declaration, element_tensor(num_entries, j)))

        # Iterate over entries in element tensor
        for (lhs, rhs) in optimized_contraction:

            # Sanity check
            ltype, i = lhs
            if ltype != 0:
                error("Expecting element tensor entry from FErari but got something else.")

            # Create name of entry
            name = element_tensor(i, j)

            # Prepare coefficents and entries for inner product
            coefficients = []
            entries = []
            for (c, rtype, k) in rhs:

                # Set coefficient
                if rtype == 0:
                    coefficients.append(c)
                    entries.append(element_tensor(k, j))
                else:
                    a = A0.secondary_multi_index.indices[k]
                    gk = geometry_tensor(j, a)
                    g_set.add(gk)
                    coefficients.append(c)
                    entries.append(gk)

            # Generate code for inner product
            value = inner_product(coefficients, entries)

            # Add declaration
            lines.append(assign(name, value))

        # Add some space
        if len(terms) > 1:
            lines.append("")

    # Sum contributions if we have more than one term
    if len(terms) > 1:
        add = format["add"]
        A0, GK, optimized_contraction = terms[0]
        for i in range(len(A0.primary_multi_index.indices)):
            name = format["element tensor"](i)
            value = add([element_tensor(i, j) for j in range(len(terms))])
            lines.append(assign(name, value))

    return "\n".join(lines)

def _generate_geometry_tensors(terms, j_set, g_set, tdim, gdim):
    "Generate code for computation of geometry tensors."

    # Prefetch formats to speed up code generation
    format_add             = format["addition"]
    format_geometry_tensor = format["geometry tensor"]
    format_scale_factor    = format["scale factor"]
    format_declaration     = format["const float declaration"]

    # Iterate over all terms
    lines = []
    offset = 0
    det_used = False

    for (i, term) in enumerate(terms):

        # Get secondary indices
        A0, GK, optimized_contraction = term
        secondary_indices = GK.secondary_multi_index.indices

        # Hack to keep old code generation based on factorization of GK
        # in case we want to reimplement factorization
        GKs = [GK]

        # Iterate over secondary indices
        for a in secondary_indices:

            # Skip code generation if term is not used
            if not format["geometry tensor"](i, a) in g_set: continue

            # Compute factorized values
            values = [_generate_entry(GK, a, offset + j, j_set, tdim, gdim) \
                          for (j, GK) in enumerate(GKs)]

            # Sum factorized values
            name = format_geometry_tensor(i, a)
            value = format_add(values)

            # Multiply with determinant factor
            det = GK.determinant
            value = _multiply_value_by_det(value, GK.determinant, len(values) > 1, j_set)
            det_used = True

            # Add code
            lines.append(format_declaration(name, value))

        # Add to offset
        offset += len(GKs)

    # Add scale factor
    if det_used:
        j_set.add(format_scale_factor) # meg says: If all values vanish, det is not used.

    return "\n".join(lines)

def _generate_entry(GK, a, i, j_set, tdim, gdim):
    "Generate code for the value of a GK entry."

    # Prefetch formats to speed up code generation
    grouping = format["grouping"]
    add      = format["addition"]
    multiply = format["multiply"]

    # Compute product of factors outside sum
    factors = _extract_factors(GK, a, None, j_set, tdim, gdim,
                              MonomialIndex.SECONDARY)

    # Compute sum of products of factors inside sum
    terms = [multiply(_extract_factors(GK, a, b, j_set, tdim, gdim,
                                       MonomialIndex.EXTERNAL))
             for b in GK.external_multi_index.indices]

    # Compute product
    if factors:
        entry = multiply(factors + [grouping(add(terms))])
    else:
        entry = add(terms)

    return entry

def _multiply_value_by_det(value, det, is_sum, j_set):
    "Generate code for multiplication of value by determinant."
    if not det.power == 0:
        J = format["det(J)"](det.restriction)
        d = [format["power"](J, det.power)]
        j_set.add(J)
    else:
        d = []
    if value == "1.0":
        v = []
    elif is_sum:
        v = [format["grouping"](value)]
    else:
        v = [value]
    return format["multiply"](d + [format["scale factor"]] + v)

def _extract_factors(GK, a, b, j_set, tdim, gdim, index_type):
    "Extract factors of given index type in GK entry."

    # Prefetch formats to speed up code generation
    coefficient = format["coefficient"]
    transform = format["transform"]

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
            # FIXME: Dimensions of J and K are transposed, what is the right thing to fix this hack?
            if t.transform_type == "J": #MonomialTransform.J:
                dim0, dim1 = gdim, tdim
            elif t.transform_type == "JINV": #MonomialTransform.JINV:
                dim0, dim1 = tdim, gdim
            else:
                error("Unknown transform type, fix this hack.")

            factors.append(transform(t.transform_type,
                                     t.index0(secondary=a, external=b),
                                     t.index1(secondary=a, external=b),
                                     dim0, dim1,
                                     t.restriction))
            j_set.add(factors[-1])

    return factors
