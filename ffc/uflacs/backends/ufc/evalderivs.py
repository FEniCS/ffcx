# -*- coding: utf-8 -*-
"""Work in progress translation of FFC evaluatebasis code to uflacs CNodes format."""

from six import string_types
import math
import numpy

from ffc.log import error
from ffc.uflacs.backends.ufc.utils import generate_error

from ffc.uflacs.backends.ufc.evaluatebasis import generate_expansion_coefficients, generate_compute_basisvalues


# Used for various indices and arrays in this file
index_type = "std::size_t"

'''
FIXME: Implement derivative mapping.
FIXME: Implement piola mapping.

// Example transformation code from evaluate_basis_derivatives:
    double transform[2][2];
    for (unsigned int j = 0; j < num_derivatives; j++)
      for (unsigned int k = 0; k < num_derivatives; k++)
        transform[j][k] = 1;
    for (unsigned int row = 0; row < num_derivatives; row++)
      for (unsigned int col = 0; col < num_derivatives; col++)
        for (unsigned int k = 0; k < n; k++)
          transform[row][col] *= Jinv[combinations[col][k]][combinations[row][k]];

// Simplified:
    double transform[max_derivatives][max_derivatives];
    for (unsigned int row = 0; row < num_derivatives; row++)
    {
      for (unsigned int col = 0; col < num_derivatives; col++)
      {
        transform[row][col] = 1.0;
        for (unsigned int k = 0; k < n; k++)
          transform[row][col] *= Jinv[combinations[col][k]][combinations[row][k]];
      }
    }

// Abstracted:
    // K = Jinv
    // C = combinations
    // Mij = product_k  K[C[i][k]][C[j][k]]
    // transform = M


// Then:
    // Reset values. Assuming that values is always an array.
    for (unsigned int r = 0; r < num_derivatives; r++)
      values[r] = 0.0;
...
      // Transform derivatives back to physical element
      for (unsigned int r = 0; r < num_derivatives; r++)
        for (unsigned int s = 0; s < num_derivatives; s++)
          values[r] += transform[r][s] * derivatives[s];
'''

def generate_evaluate_reference_basis_derivatives(L, data, parameters):
    # Cutoff for feature to disable generation of this code (consider removing after benchmarking final result)
    if isinstance(data, string_types):
        msg = "evaluate_reference_basis_derivatives: %s" % (data,)
        return generate_error(L, msg, parameters["convert_exceptions_to_warnings"])

    # Get some known dimensions
    element_cellname = data["cellname"]
    tdim = data["topological_dimension"]
    max_degree = data["max_degree"]
    reference_value_size = data["reference_value_size"]
    num_dofs = len(data["dofs_data"])

    # Output argument
    reference_values = L.Symbol("reference_values")

    # Input arguments
    order = L.Symbol("order")
    num_points = L.Symbol("num_points")
    X = L.Symbol("X")

    # Loop indices
    ip = L.Symbol("ip") # point
    idof = L.Symbol("i")   # dof
    c = L.Symbol("c")   # component
    r = L.Symbol("r")   # derivative number

    # Define symbol for number of derivatives of given order
    num_derivatives = L.Symbol("num_derivatives")
    reference_values_size = num_points * num_dofs * num_derivatives * reference_value_size

    # FIXME: validate these dimensions
    ref_values = L.FlattenedArray(reference_values,
        dims=(num_points, num_dofs, num_derivatives, reference_value_size))
    # From evaluatebasis.py:
    #ref_values = L.FlattenedArray(reference_values, dims=(num_points, num_dofs, reference_value_size))

    # Initialization (zeroing) and cutoffs outside valid range of orders
    setup_code = [
        # Cutoff to evaluate_basis for order 0
        L.If(L.EQ(order, 0), [
            L.Call("evaluate_reference_basis", (reference_values, num_points, X)),
            L.Return()
            ]),
        # Compute number of derivatives of this order
        L.VariableDecl("const " + index_type, num_derivatives,
                       value=L.Call("std::pow", (tdim, order))),
        # Reset output values to zero
        L.MemZero(reference_values, reference_values_size),
        # Cutoff for higher order than we have
        L.If(L.GT(order, max_degree),
             L.Return()),
        ]

    # If max_degree is zero, we don't need to generate any more code
    if max_degree == 0:
        return setup_code

    # Tabulate dmats tables for all dofs and all derivative directions
    dmats_names, dmats_code = generate_tabulate_dmats(L, data["dofs_data"])

    # Generate code with static tables of expansion coefficients
    tables_code, coefficients_for_dof = generate_expansion_coefficients(L, data["dofs_data"])


    # Generate code to compute tables of basisvalues
    basisvalues_code, basisvalues_for_degree, need_fiat_coordinates = \
        generate_compute_basisvalues(L, data["dofs_data"], element_cellname, tdim, X, ip)

    # Accumulate products of basisvalues and coefficients into values
    accumulation_code = generate_accumulation_code(L,
        data["dofs_data"],
        num_derivatives,
        coefficients_for_dof)

    # Generate geo code.
    geometry_code = [] # _geometry_related_code(L, data, tdim, gdim, element_cellname)  # FIXME

    # Generate all possible combinations of derivatives.
    combinations_code, combinations = _generate_combinations(L, tdim, max_degree, order, num_derivatives)

    # Define symbols for variables populated inside dof switch
    derivatives = L.Symbol("derivatives")
    reference_offset = L.Symbol("reference_offset")
    num_components = L.Symbol("num_components")

    # Get number of components of each basis function (>1 for dofs of piola mapped subelements)
    num_components_values = [dof_data["num_components"] for dof_data in data["dofs_data"]]

    # Offset into parent mixed element to first component of each basis function
    reference_offset_values = [dof_data["reference_offset"] for dof_data in data["dofs_data"]]

    # Max dimensions for the reference derivatives for each dof
    max_num_derivatives = tdim**max_degree
    max_num_components = max(num_components_values)

    # Add constant tables of these numbers
    tables_code += [
        L.ArrayDecl("const " + index_type, reference_offset, num_dofs, values=reference_offset_values),
        L.ArrayDecl("const " + index_type, num_components, num_dofs, values=num_components_values),
    ]

    # Access reference derivatives compactly
    derivs = L.FlattenedArray(derivatives, dims=(num_components[idof], num_derivatives))

    # Compute aux[i] = sum_j dmats[i][j] * basisvalues[j] for all unique combinations
    all_aux_computation = []
    aux_names = {}
    aux_for_dof = {}
    for i_dof, dof_data in enumerate(data["dofs_data"]):
        embedded_degree = dof_data["embedded_degree"]
        basisvalues = basisvalues_for_degree[embedded_degree]

        key = (dmats_names[i_dof].name, basisvalues.name)
        aux = aux_names.get(key)
        if aux is None:
            aux_computation, aux = _compute_aux_dmats_basisvalues_products(L, dof_data,
                i_dof, order, num_derivatives, combinations,
                dmats_names, basisvalues)
            aux_names[key] = aux
            all_aux_computation += aux_computation

        aux_for_dof[i_dof] = aux

    # Create code for all basis values (dofs).
    dof_cases = []
    for i_dof, dof_data in enumerate(data["dofs_data"]):
        embedded_degree = dof_data["embedded_degree"]
        basisvalues = basisvalues_for_degree[embedded_degree]

        # Compute the derivatives of the basisfunctions on the reference (FIAT) element,
        # as the dot product of the new coefficients and basisvalues.
        case_code = _compute_reference_derivatives(L, dof_data,
                        i_dof, num_derivatives, derivs,
                        coefficients_for_dof, aux_for_dof)

        dof_cases.append((i_dof, case_code))

    # Loop over all dofs, entering a different switch case in each iteration.
    # This is a legacy from the previous implementation where the loop
    # was in a separate function and all the setup above was also repeated
    # in a call for each dof.
    # TODO: Further refactoring is needed to improve on this situation,
    # but at least it's better than before. There's probably more code and
    # computations that can be shared between dofs, and this would probably
    # be easier to fix if mixed elements were handled separately!
    dof_loop_code = [
        L.Comment("Loop over all dofs"),
        L.ForRange(idof, 0, num_dofs, index_type=index_type, body=[
            L.ArrayDecl("double", derivatives, max_num_components * max_num_derivatives),
            L.Switch(idof, dof_cases),
            L.ForRange(r, 0, num_derivatives, index_type=index_type, body=[
                L.ForRange(c, 0, num_components[idof], index_type=index_type, body=[
                    L.Assign(ref_values[ip][idof][r][reference_offset[idof] + c], derivs[c][r]),  # FIXME: validate ref_values dims
                ]),
            ])
        ]),
    ]

    # Define loop over points
    final_loop_code = [
        L.ForRange(ip, 0, num_points, index_type=index_type, body=
            basisvalues_code
            + all_aux_computation
            + dof_loop_code
            )
        ]

    # Stitch it all together
    code = (
        setup_code
        + dmats_code
        + tables_code
        + geometry_code
        + combinations_code
#        + dof_loop_code
        + final_loop_code
        )
    return code


def generate_tabulate_dmats(L, dofs_data):
    "Tabulate the derivatives of the polynomial base"

    alignas = 32

    # Emit code for the dmats we've actually used
    dmats_code = [L.Comment("Tables of derivatives of the polynomial base (transpose).")]

    dmats_names = []

    all_matrices = []

    for idof, dof_data in enumerate(dofs_data):
        # Get derivative matrices (coefficients) of basis functions, computed by FIAT at compile time.
        derivative_matrices = dof_data["dmats"]
        num_mats = len(derivative_matrices)
        num_members = dof_data["num_expansion_members"]

        # Generate tables for each spatial direction.
        matrix = numpy.zeros((num_mats, num_members, num_members))
        for i, dmat in enumerate(derivative_matrices):
            # Extract derivatives for current direction
            # (take transpose, FIAT_NEW PolynomialSet.tabulate()).
            matrix[i,...] = numpy.transpose(dmat)

        # O(n^2) matrix matching...
        name = None
        for oldname, oldmatrix in all_matrices:
            if matrix.shape == oldmatrix.shape and numpy.allclose(matrix, oldmatrix):
                name = oldname
                break

        if name is None:
            # Define variable name for coefficients for this dof
            name = L.Symbol("dmats%d" % (idof,))
            all_matrices.append((name, matrix))

            # Declare new dmats table with unique values
            decl = L.ArrayDecl("static const double", name, (num_mats, num_members, num_members),
                               values=matrix, alignas=alignas)
            dmats_code.append(decl)

        # Append name for each dof
        dmats_names.append(name)

    return dmats_names, dmats_code


def _generate_combinations(L, tdim, max_degree, order, num_derivatives):
    max_num_derivatives = tdim**max_degree
    combinations = L.Symbol("combinations")

    # This precomputes the combinations for each order and stores in code as table
    # Python equivalent precomputed for each valid order:
    combinations_shape = (max_degree, max_num_derivatives, max_degree)
    all_combinations = numpy.zeros(combinations_shape, dtype=int)
    for q in range(1, max_degree+1):
        for row in range(1, max_num_derivatives):
            for num in range(0, row):
                for col in range(q-1, -1, -1):
                    if all_combinations[q-1][row][col] > tdim - 2:
                        all_combinations[q-1][row][col] = 0
                    else:
                        all_combinations[q-1][row][col] += 1
                        break
    code = [
        L.Comment("Precomputed combinations"),
        L.ArrayDecl("const " + index_type, combinations, combinations_shape, values=all_combinations),
        ]
    # Select the right order for further access
    combinations = combinations[order-1]

    return code, combinations


def generate_accumulation_code(L, dofs_data, num_derivatives, coefficients_for_dof):
    code = []
    # FIXME
    return code


def _compute_aux_dmats_basisvalues_products(L, dof_data,
        idof, order, num_derivatives, combinations,
        dmats_names, basisvalues):
    """Deprecated comment after refactoring, but contents are still relevant:

    Compute derivatives on the reference element by recursively multiply coefficients with
    the relevant derivatives of the polynomial base until the requested order of derivatives
    has been reached. After this take the dot product with the basisvalues."""

    r = L.Symbol("r")
    s = L.Symbol("s")
    t = L.Symbol("t")
    u = L.Symbol("u")
    tu = L.Symbol("tu")

    # Get number of components.
    num_components = dof_data["num_components"]  # FIXME: Is this the number of components of the subelement?

    # Get shape of derivative matrix (they should all have the same shape) and
    # verify that it is a square matrix.
    shape_dmats = numpy.shape(dof_data["dmats"][0])
    if shape_dmats[0] != shape_dmats[1]:
        error("Something is wrong with the dmats:\n%s" % str(dof_data["dmats"]))

    # Declare matrix of dmats (which will hold the matrix product of all combinations)
    # and dmats_old which is needed in order to perform the matrix product.
    dmats = L.Symbol("dmats")
    temp_dmats_declarations = [
        L.Comment("Declare derivative matrix (of polynomial basis)."),
        L.ArrayDecl("double", dmats, shape_dmats),
        ]

    # Compute dmats as a recursive matrix product
    dmats_computation = _compute_dmats(L, len(dof_data["dmats"]), shape_dmats,
                                       dmats, dmats_names, idof, order, combinations,
                                       [s, t, u, tu], r)

    # Accumulate aux = dmats * basisvalues
    aux = L.Symbol("aux_%s_%s" % (dmats_names[idof], basisvalues))
    aux_declaration = L.ArrayDecl("double", aux, shape_dmats[0], values=0)
    aux_accumulation = [
        L.ForRange(s, 0, shape_dmats[0], index_type=index_type, body=[
            L.ForRange(t, 0, shape_dmats[1], index_type=index_type, body=[
                L.AssignAdd(aux[s], dmats[s, t] * basisvalues[t])
            ])
        ])
    ]

    # Generate loop over number of derivatives.
    # Loop all derivatives and compute value of the derivative as:
    # deriv_on_ref[r] = coeff[dof][s]*dmat[s][t]*basis[t]
    aux_computation = [
        aux_declaration,
        L.ForRange(r, 0, num_derivatives, index_type=index_type, body=
            temp_dmats_declarations
            + dmats_computation
            + aux_accumulation
        )
    ]
    return aux_computation, aux


def _compute_reference_derivatives(L, dof_data, idof,
        num_derivatives, derivatives,
        coefficients_for_dof, aux_for_dof):
    """Deprecated comment after refactoring, but contents are still relevant:

    Compute derivatives on the reference element by recursively multiply coefficients with
    the relevant derivatives of the polynomial base until the requested order of derivatives
    has been reached. After this take the dot product with the basisvalues."""

    r = L.Symbol("r")
    s = L.Symbol("s")

    # Get shape of derivative matrix (they should all have the same shape)
    shape_dmats = numpy.shape(dof_data["dmats"][0])

    # Get number of components of this basis function
    # (>1 for dofs of piola mapped subelements)
    num_components = dof_data["num_components"]

    # Compute derivatives for all derivatives and components for this particular idof
    code = [
        L.Comment("Compute reference derivatives for dof %d." % idof),
        # Accumulate sum_s coefficients[s] * aux[s]
        L.ForRange(r, 0, num_derivatives, index_type=index_type, body=[
            # Unrolled loop over components of basis function
            L.Assign(derivatives[c][r], 0.0)
            for c in range(num_components)
            ] + [
            L.ForRange(s, 0, shape_dmats[0], index_type=index_type, body=[
                # Unrolled loop over components of basis function
                L.AssignAdd(derivatives[c][r],
                            coefficients_for_dof[idof][c][s] * aux_for_dof[idof][s])
                for c in range(num_components)
            ])
        ])
    ]
    return code


def _compute_dmats(L, num_dmats, shape_dmats,
                   dmats, dmats_names, idof, order, combinations,
                   available_indices, deriv_index):
    "Compute values of dmats as a matrix product."
    s, t, u, tu = available_indices

    dmats_old = L.Symbol("dmats_old")

    # Create dmats matrix by multiplication
    comb = L.Symbol("comb")

    dmats_name = dmats_names[idof]

    # Local helper function for loops over t,u < shape_dmats
    def tu_loops(body):
        return L.ForRange(t, 0, shape_dmats[0], index_type=index_type, body=
                   L.ForRange(u, 0, shape_dmats[1], index_type=index_type, body=
                       body
                   )
               )

    code = [
        L.Comment("Initialize dmats."),
        L.VariableDecl(index_type, comb, combinations[deriv_index, 0]),
        tu_loops(
            L.Assign(dmats[t, u], dmats_name[comb, t, u])
        ),
        L.Comment("Looping derivative order to generate dmats."),
        L.ForRange(s, 1, order, index_type=index_type, body=[
            L.Comment("Store previous dmats matrix."),
            L.ArrayDecl("double", dmats_old, shape_dmats),
            tu_loops(
                L.Assign(dmats_old[t, u], dmats[t, u])
            ),
            L.Comment("Resetting dmats."),
            tu_loops(
                L.Assign(dmats[t, u], L.LiteralFloat(0.0))
            ),
            L.Comment("Update dmats using an inner product."),
            L.Assign(comb, combinations[deriv_index, s]),
            tu_loops(
                L.ForRange(tu, 0, shape_dmats[0], index_type=index_type, body=
                        L.AssignAdd(dmats[t, u], dmats_name[comb, t, tu] * dmats_old[tu, u])
                )
            ),
        ])
    ]
    return code

'''

                    // Actually instead want to do this for each (unique dmats%d, basisvalues%d) combo
                    double tmp[num_derivatives][num_members] = {};
                    for (std::size_t r = 0; r < num_derivatives; ++r)
                    {
                        // Declare derivative matrix (of polynomial basis).
                        double dmats[3][3] =
                            { { 1.0, 0.0, 0.0 },
                              { 0.0, 1.0, 0.0 },
                              { 0.0, 0.0, 1.0 } };
                        // Looping derivative order to generate dmats.
                        for (std::size_t s = 0; s < order; ++s)
                        {
                            // Store previous dmats matrix.
                            double dmats_old[3][3];
                            for (std::size_t t = 0; t < 3; ++t)
                                for (std::size_t u = 0; u < 3; ++u)
                                    dmats_old[t][u] = dmats[t][u];
                            // Resetting dmats.
                            for (std::size_t t = 0; t < 3; ++t)
                                for (std::size_t u = 0; u < 3; ++u)
                                    dmats[t][u] = 0.0;
                            // Update dmats using an inner product.
                            for (std::size_t t = 0; t < 3; ++t)
                                for (std::size_t u = 0; u < 3; ++u)
                                    for (std::size_t tu = 0; tu < 3; ++tu)
                                        dmats[t][u] += dmats_A[combinations[order - 1][r][s]][t][tu] * dmats_old[tu][u];
                        }

                        for (std::size_t s = 0; s < 3; ++s)
                            for (std::size_t t = 0; t < 3; ++t)
                                tmp[r][s] += dmats[s][t] * basisvalues1[t];
                    }


'''
'''
        // Loop over all dofs
        for (std::size_t i = 0; i < 3; ++i)
        {
            double derivatives[2] = {};
            switch (i)
            {
            case 0:
                // Only need this inside the per-dof code:
                for (std::size_t r = 0; r < num_derivatives; ++r)
                  for (std::size_t s = 0; s < 3; ++s)
                    derivatives[r] += coefficients0[0][s] * tmp[r][s];
                break;
'''
