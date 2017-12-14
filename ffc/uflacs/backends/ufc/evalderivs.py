# -*- coding: utf-8 -*-
"""Work in progress translation of FFC evaluatebasis code to uflacs CNodes format."""

import numpy

from ffc.log import error
from ffc.uflacs.backends.ufc.utils import generate_error

from ffc.uflacs.backends.ufc.evaluatebasis import generate_expansion_coefficients, generate_compute_basisvalues


# Used for various indices and arrays in this file
index_type = "std::size_t"


def generate_evaluate_reference_basis_derivatives(L, data, parameters):
    # Cutoff for feature to disable generation of this code (consider removing after benchmarking final result)
    if isinstance(data, str):
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
                                  dims=(num_points, num_dofs, num_derivatives,
                                        reference_value_size))
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

    # Create code for all basis values (dofs).
    dof_cases = []
    for i_dof, dof_data in enumerate(data["dofs_data"]):

        embedded_degree = dof_data["embedded_degree"]
        basisvalues = basisvalues_for_degree[embedded_degree]

        shape_dmats = numpy.shape(dof_data["dmats"][0])
        if shape_dmats[0] != shape_dmats[1]:
            error("Something is wrong with the dmats:\n%s" % str(dof_data["dmats"]))

        aux = L.Symbol("aux")
        dmats = L.Symbol("dmats")
        dmats_old = L.Symbol("dmats_old")
        dmats_name = dmats_names[i_dof]

        # Create dmats matrix by multiplication
        comb = L.Symbol("comb")
        s = L.Symbol("s")
        t = L.Symbol("t")
        u = L.Symbol("u")
        tu = L.Symbol("tu")
        aux_computation_code = [ L.ArrayDecl("double", aux, shape_dmats[0], values=0),
                            L.Comment("Declare derivative matrix (of polynomial basis)."),
                            L.ArrayDecl("double", dmats, shape_dmats, values=0),
                            L.Comment("Initialize dmats."),
                            L.VariableDecl(index_type, comb, combinations[r, 0]),
                            L.MemCopy(L.AddressOf(dmats_name[comb][0][0]),  L.AddressOf(dmats[0][0]), shape_dmats[0]*shape_dmats[1]),
                            L.Comment("Looping derivative order to generate dmats."),
                            L.ForRange(s, 1, order, index_type=index_type, body=[
                                L.Comment("Store previous dmats matrix."),
                                L.ArrayDecl("double", dmats_old, shape_dmats),
                                L.MemCopy(L.AddressOf(dmats[0][0]),  L.AddressOf(dmats_old[0][0]), shape_dmats[0]*shape_dmats[1]),
                                L.Comment("Resetting dmats."),
                                L.MemZero(L.AddressOf(dmats[0][0]), shape_dmats[0]*shape_dmats[1]),
                                L.Comment("Update dmats using an inner product."),
                                L.Assign(comb, combinations[r, s]),
                                L.ForRange(t, 0, shape_dmats[0], index_type=index_type, body=
                                   L.ForRange(u, 0, shape_dmats[1], index_type=index_type, body=
                                      L.ForRange(tu, 0, shape_dmats[0], index_type=index_type, body=
                                         L.AssignAdd(dmats[t, u], dmats_name[comb, t, tu] * dmats_old[tu, u]))))]),
                            L.ForRange(s, 0, shape_dmats[0], index_type=index_type, body=
                               L.ForRange(t, 0, shape_dmats[1], index_type=index_type, body=
                                          L.AssignAdd(aux[s], dmats[s, t] * basisvalues[t])))
        ]

        # Unrolled loop over components of basis function
        n = dof_data["num_components"]
        compute_ref_derivs_code = [L.Assign(derivs[cc][r], 0.0) for cc in range(n)]

        compute_ref_derivs_code += [L.ForRange(s, 0, shape_dmats[0], index_type=index_type, body=
                                               [L.AssignAdd(derivs[cc][r], coefficients_for_dof[i_dof][cc][s] * aux[s]) for cc in range(n)])]

        embedded_degree = dof_data["embedded_degree"]
        basisvalues = basisvalues_for_degree[embedded_degree]

        # Compute the derivatives of the basisfunctions on the reference (FIAT) element,
        # as the dot product of the new coefficients and basisvalues.

        case_code =  [L.Comment("Compute reference derivatives for dof %d." % i_dof),
                      # Accumulate sum_s coefficients[s] * aux[s]
                      L.ForRange(r, 0, num_derivatives, index_type=index_type, body=[
                          aux_computation_code,
                          compute_ref_derivs_code
                      ])]

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
            L.ArrayDecl("double", derivatives, max_num_components * max_num_derivatives, 0.0),
            L.Switch(idof, dof_cases),
            L.ForRange(r, 0, num_derivatives, index_type=index_type, body= [
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
            + dof_loop_code
            )
        ]

    # Stitch it all together
    code = (
        setup_code
        + dmats_code
        + tables_code
        + combinations_code
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

        # TODO: Use precision from parameters here
        from ffc.uflacs.elementtables import clamp_table_small_numbers
        matrix = clamp_table_small_numbers(matrix)

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


def _generate_combinations(L, tdim, max_degree, order, num_derivatives, suffix=""):
    max_num_derivatives = tdim**max_degree
    combinations = L.Symbol("combinations" + suffix)

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
