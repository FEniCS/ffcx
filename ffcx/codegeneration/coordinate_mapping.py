# Copyright (C) 2015-2017 Martin Sandve Alnæs
#
# This file is part of FFCX.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import ffcx.codegeneration.coordinate_mapping_template as ufc_coordinate_mapping

# TODO: Test everything here! Cover all combinations of gdim,tdim=1,2,3!

index_type = "int64_t"

# Code generation utilities:


def generate_compute_ATA(L, ATA, A, m, n, index_prefix=""):
    """Generate code to declare and compute ATA[i,j] = sum_k A[k,i]*A[k,j] with given A shaped (m,n)."""
    # Loop indices
    i = L.Symbol(index_prefix + "i")
    j = L.Symbol(index_prefix + "j")
    k = L.Symbol(index_prefix + "k")

    # Build A^T*A matrix
    code = [
        L.ArrayDecl("double", ATA, (n, n), values=0),
        L.ForRanges(
            (i, 0, n), (j, 0, n), (k, 0, m),
            index_type=index_type,
            body=L.AssignAdd(ATA[i, j], A[k, i] * A[k, j])),
    ]
    return L.StatementList(code)


# Inline math expressions:


def det_22(B, i, j, k, l):
    return B[i, k] * B[j, l] - B[i, l] * B[j, k]


def codet_nn(A, rows, cols):
    n = len(rows)
    if n == 2:
        return det_22(A, rows[0], rows[1], cols[0], cols[1])
    else:
        r = rows[0]
        subrows = rows[1:]
        parts = []
        for i, c in enumerate(cols):
            subcols = cols[i + 1:] + cols[:i]
            parts.append(A[r, c] * codet_nn(A, subrows, subcols))
        return sum(parts[1:], parts[0])


def det_nn(A, n):
    if n == 1:
        return A[0, 0]
    else:
        ns = list(range(n))
        return codet_nn(A, ns, ns)


def pdet_m1(L, A, m):
    # Special inlined case 1xm for simpler expression
    A2 = A[0, 0] * A[0, 0]
    for i in range(1, m):
        A2 = A2 + A[i, 0] * A[i, 0]
    return L.Sqrt(A2)


def adj_expr_2x2(A):
    return [[A[1, 1], -A[0, 1]], [-A[1, 0], A[0, 0]]]


def adj_expr_3x3(A):
    return [[
        A[2, 2] * A[1, 1] - A[1, 2] * A[2, 1], A[0, 2] * A[2, 1] - A[0, 1] * A[2, 2],
        A[0, 1] * A[1, 2] - A[0, 2] * A[1, 1]
    ], [
        A[1, 2] * A[2, 0] - A[2, 2] * A[1, 0], A[2, 2] * A[0, 0] - A[0, 2] * A[2, 0],
        A[0, 2] * A[1, 0] - A[1, 2] * A[0, 0]
    ], [
        A[1, 0] * A[2, 1] - A[2, 0] * A[1, 1], A[0, 1] * A[2, 0] - A[0, 0] * A[2, 1],
        A[0, 0] * A[1, 1] - A[0, 1] * A[1, 0]
    ]]


def generate_assign_inverse(L, K, J, detJ, gdim, tdim):
    if gdim == tdim:
        if gdim == 1:
            return L.Assign(K[0, 0], L.LiteralFloat(1.0) / J[0, 0])
        elif gdim == 2:
            adj_values = adj_expr_2x2(J)
        elif gdim == 3:
            adj_values = adj_expr_3x3(J)
        else:
            raise RuntimeError("Not implemented.")
        return L.StatementList(
            [L.Assign(K[j, i], adj_values[j][i] / detJ) for j in range(tdim) for i in range(gdim)])
    else:
        if tdim == 1:
            # Simpler special case for embedded 1d
            prods = [J[i, 0] * J[i, 0] for i in range(gdim)]
            s = sum(prods[1:], prods[0])
            return L.StatementList([L.Assign(K[0, i], J[i, 0] / s) for i in range(gdim)])
        else:
            # Generic formulation of Penrose-Moore pseudo-inverse of J: (J.T*J)^-1 * J.T
            i = L.Symbol("i")
            j = L.Symbol("j")
            k = L.Symbol("k")
            JTJ = L.Symbol("JTJ")
            JTJf = L.FlattenedArray(JTJ, dims=(tdim, tdim))
            detJTJ = detJ * detJ
            JTJinv = L.Symbol("JTJinv")
            JTJinvf = L.FlattenedArray(JTJinv, dims=(tdim, tdim))
            code = [
                L.Comment("Compute J^T J"),
                L.ArrayDecl("double", JTJ, (tdim * tdim, ), values=0),
                L.ForRanges(
                    (k, 0, tdim), (j, 0, tdim), (i, 0, gdim),
                    index_type=index_type,
                    body=L.AssignAdd(JTJf[k, j], J[i, k] * J[i, j])),
                L.Comment("Compute inverse(J^T J)"),
                L.ArrayDecl("double", JTJinv, (tdim * tdim, )),
                generate_assign_inverse(L, JTJinvf, JTJf, detJTJ, tdim, tdim),
                L.Comment("Compute K = inverse(J^T J) * J"),
                L.ForRanges(
                    (k, 0, tdim), (i, 0, gdim),
                    index_type=index_type,
                    body=L.Assign(K[k, i], L.LiteralFloat(0.0))),
                L.ForRanges(
                    (k, 0, tdim), (i, 0, gdim), (j, 0, tdim),
                    index_type=index_type,
                    body=L.AssignAdd(K[k, i], JTJinvf[k, j] * J[i, j])),
            ]
            return L.StatementList(code)


def evaluate_reference_basis_declaration(L, ir):
    scalar_coordinate_element_classname = ir.scalar_coordinate_finite_element_classname
    code = """
int evaluate_reference_basis_{}(double* restrict reference_values,
    int num_points, const double* restrict X);
""".format(scalar_coordinate_element_classname)
    return code


def compute_physical_coordinates(L, ir):
    num_dofs = ir.num_scalar_coordinate_element_dofs
    scalar_coordinate_element_classname = ir.scalar_coordinate_finite_element_classname

    # Dimensions
    gdim = ir.geometric_dimension
    tdim = ir.topological_dimension
    num_points = L.Symbol("num_points")

    # Loop indices
    ip = L.Symbol("ip")
    i = L.Symbol("i")
    d = L.Symbol("d")

    # Input cell data
    coordinate_dofs = L.FlattenedArray(L.Symbol("coordinate_dofs"), dims=(num_dofs, gdim))

    # Output geometry
    x = L.FlattenedArray(L.Symbol("x"), dims=(num_points, gdim))

    # Input geometry
    X = L.FlattenedArray(L.Symbol("X"), dims=(num_points, tdim))

    # Computing table one point at a time instead of
    # using num_points to avoid dynamic allocation
    one_point = 1

    # Symbols for local basis values table
    phi_sym = L.Symbol("phi")

    # NB! Must match array layout of evaluate_reference_basis
    phi = L.FlattenedArray(phi_sym, dims=(num_dofs, ))

    # For each point, compute basis values and accumulate into the right x
    code = [
        L.ArrayDecl("double", phi_sym, (one_point * num_dofs, )),
        L.ForRange(i, 0, num_points * gdim, index_type=index_type, body=L.Assign(x.array[i], 0.0)),
        L.ForRange(
            ip,
            0,
            num_points,
            index_type=index_type,
            body=[
                L.Comment("Compute basis values of coordinate element"),
                L.Call("evaluate_reference_basis_{}".format(scalar_coordinate_element_classname),
                       (phi_sym, 1, L.AddressOf(X[ip, 0]))),
                L.Comment("Compute x"),
                L.ForRanges(
                    (i, 0, gdim), (d, 0, num_dofs),
                    index_type=index_type,
                    body=L.AssignAdd(x[ip, i], coordinate_dofs[d, i] * phi[d]))
            ]),
    ]

    return code


def compute_reference_geometry(L, ir):
    degree = ir.coordinate_element_degree
    if degree == 1:
        # Special case optimized for affine mesh (possibly room for
        # further optimization)
        return _compute_reference_coordinates_affine(L, ir, output_all=True)
    else:
        # General case with newton loop to solve F(X) = x(X) - x0 = 0
        return _compute_reference_coordinates_newton(L, ir, output_all=True)


# TODO: Maybe we don't need this version, see what we need in dolfinx first
def compute_reference_coordinates(L, ir):
    degree = ir.coordinate_element_degree
    if degree == 1:
        # Special case optimized for affine mesh (possibly room for
        # further optimization)
        return _compute_reference_coordinates_affine(L, ir)
    else:
        # General case with newton loop to solve F(X) = x(X) - x0 = 0
        return _compute_reference_coordinates_newton(L, ir)


def _compute_reference_coordinates_affine(L, ir, output_all=False):
    # Class name
    classname = ir.name

    # Dimensions
    gdim = ir.geometric_dimension
    tdim = ir.topological_dimension
    num_points = L.Symbol("num_points")

    # Number of dofs for a scalar component
    num_dofs = ir.num_scalar_coordinate_element_dofs

    # Loop indices
    ip = L.Symbol("ip")  # point
    i = L.Symbol("i")  # gdim
    j = L.Symbol("j")  # tdim
    k = L.Symbol("k")  # sum iteration
    iz = L.Symbol("l")  # zero X array

    # Output geometry
    X = L.FlattenedArray(L.Symbol("X"), dims=(num_points, tdim))
    # if output_all, this is also output:
    Jsym = L.Symbol("J")
    detJsym = L.Symbol("detJ")
    Ksym = L.Symbol("K")

    # Input geometry
    x = L.FlattenedArray(L.Symbol("x"), dims=(num_points, gdim))

    # Input cell data
    coordinate_dofs = L.FlattenedArray(L.Symbol("coordinate_dofs"), dims=(num_dofs, gdim))

    init_input = [
        L.ForRange(
            iz, 0, num_points * tdim, index_type=index_type, body=L.Assign(X.array[iz], 0.0))
    ]

    if output_all:
        decls = []
    else:
        decls = [
            L.ArrayDecl("double", Jsym, sizes=(gdim * tdim, )),
            L.ArrayDecl("double", detJsym, sizes=(1, )),
            L.ArrayDecl("double", Ksym, sizes=(tdim * gdim, )),
        ]

    # Tables of coordinate basis function values and derivatives at X=0
    # and X=midpoint available through ir. This is useful in several
    # geometry functions.
    tables = ir.tables

    # Check the table shapes against our expectations
    x_table = tables["x0"]
    J_table = tables["J0"]
    assert x_table.shape == (num_dofs, )
    assert J_table.shape == (tdim, num_dofs)

    # TODO: Use epsilon parameter here?
    # TODO: Move to a more 'neutral' utility file
    from ffcx.ir.uflacs.elementtables import clamp_table_small_numbers
    x_table = clamp_table_small_numbers(x_table)
    J_table = clamp_table_small_numbers(J_table)

    # Table symbols
    phi_X0 = L.Symbol("phi_X0")
    dphi_X0 = L.Symbol("dphi_X0")

    # Table declarations
    table_decls = [
        L.ArrayDecl("const double", phi_X0, sizes=x_table.shape, values=x_table),
        L.ArrayDecl("const double", dphi_X0, sizes=J_table.shape, values=J_table),
    ]

    # Compute x0 = x(X=0) (optimized by precomputing basis at X=0)
    x0 = L.Symbol("x0")
    compute_x0 = [
        L.ArrayDecl("double", x0, sizes=(gdim, ), values=0),
        L.ForRanges(
            (i, 0, gdim), (k, 0, num_dofs),
            index_type=index_type,
            body=L.AssignAdd(x0[i], coordinate_dofs[k, i] * phi_X0[k])),
    ]

    # For more convenient indexing
    J = L.FlattenedArray(Jsym, dims=(gdim, tdim))
    K = L.FlattenedArray(Ksym, dims=(tdim, gdim))

    # Compute J = J(X=0) (optimized by precomputing basis at X=0)
    compute_J0 = [
        L.ForRanges(
            (i, 0, gdim), (j, 0, tdim),
            index_type=index_type,
            body=[
                L.Assign(J[i, j], 0.0),
                L.ForRange(
                    k,
                    0,
                    num_dofs,
                    index_type=index_type,
                    body=L.AssignAdd(J[i, j], coordinate_dofs[k, i] * dphi_X0[j, k]))
            ]),
    ]

    # Compute K = inv(J) (and intermediate value det(J))
    compute_K0 = [
        L.Call("compute_jacobian_determinants_{}".format(classname),
               (detJsym, 1, Jsym)),
        L.Call("compute_jacobian_inverses_{}".format(classname), (Ksym, 1, Jsym, detJsym)),
    ]

    # Compute X = K0*(x-x0) for each physical point x
    compute_X = [
        L.ForRanges(
            (ip, 0, num_points), (j, 0, tdim), (i, 0, gdim),
            index_type=index_type,
            body=L.AssignAdd(X[ip, j], K[j, i] * (x[ip, i] - x0[i])))
    ]

    # Stitch it together
    code = init_input + table_decls + decls + compute_x0 + compute_J0 + compute_K0 + compute_X
    return code


def _compute_reference_coordinates_newton(L, ir, output_all=False):
    """Solves x(X) = x0 for X.

    Find X such that, given x0,

        F(X) = x(X) - x0 = 0

    Newton iteration is:

        X^0 = midpoint of reference cell
        until convergence:
            dF/dX(X^k) dX^k = -F(X^k)
            X^{k+1} = X^k + dX^k

    The Jacobian is just the usual Jacobian of the geometry mapping:

        dF/dX = dx/dX = J

    and its inverse is the usual K = J^-1.

    Because J is small, we can invert it directly and rewrite

        Jk dX = -(xk - x0) = x0 - xk

    into

        Xk = midpoint of reference cell
        for k in range(maxiter):
            xk = x(Xk)
            Jk = J(Xk)
            Kk = inverse(Jk)
            dX = Kk (x0 - xk)
            Xk += dX
            if dX sufficiently small: break
    """
    # Dimensions
    gdim = ir.geometric_dimension
    tdim = ir.topological_dimension
    cellname = ir.cell_shape
    num_points = L.Symbol("num_points")

    degree = ir.coordinate_element_degree

    # Computing table one point at a time instead of vectorized over
    # num_points will allow skipping dynamic allocation
    one_point = 1

    # Loop indices
    ip = L.Symbol("ip")  # point
    i = L.Symbol("i")  # gdim
    j = L.Symbol("j")  # tdim
    k = L.Symbol("k")  # iteration

    # Input cell data
    coordinate_dofs = L.Symbol("coordinate_dofs")

    # Output geometry
    X = L.FlattenedArray(L.Symbol("X"), dims=(num_points, tdim))

    # Input geometry
    x = L.FlattenedArray(L.Symbol("x"), dims=(num_points, gdim))

    # Find X=Xk such that xk(Xk) = x or F(Xk) = xk(Xk) - x = 0.
    # Newtons method is then:
    #   X0 = cell midpoint
    #   dF/dX dX = -F
    #   Xk = Xk + dX
    # Note that
    #   dF/dX = dx/dX = J
    # giving
    #   inverse(dF/dX) = K
    # such that dX = -K*(xk(Xk) - x)

    # Newtons method is then:
    # for each ip:
    #   xgoal = x[ip]
    #   X0 = cell midpoint
    #   until convergence:
    #       K = K(Xk)
    #       dX = K*(xk(Xk) - x)
    #       Xk -= dX
    #   X[ip] = Xk

    # Symbols for arrays used below
    Xk = L.Symbol("Xk")
    xgoal = L.Symbol("xgoal")
    dX = L.Symbol("dX")
    xk = L.Symbol("xk")
    J = L.Symbol("J")
    detJ = L.Symbol("detJ")
    K = L.Symbol("K")

    xm = L.Symbol("xm")
    Km = L.Symbol("Km")

    # Symbol for ufc_geometry cell midpoint definition
    Xm = L.Symbol("%s_midpoint" % cellname)

    # Variables for stopping criteria
    # TODO: Check if these are good convergence criteria,
    #       e.g. is epsilon=1e-6 and iterations=degree sufficient?
    max_iter = L.LiteralInt(degree)
    epsilon = L.LiteralFloat(1e-6)
    # TODO: Could also easily make criteria input if desired
    # max_iter = L.Symbol("iterations")
    # epsilon = L.Symbol("epsilon")
    dX2 = L.Symbol("dX2")

    # Wrap K as flattened array for convenient indexing Kf[j,i]
    Kf = L.FlattenedArray(K, dims=(tdim, gdim))
    Kmf = L.FlattenedArray(Km, dims=(tdim, gdim))

    decls = [
        L.Comment("Declare intermediate arrays to hold results of compute_geometry call"),
        L.ArrayDecl("double", xk, (gdim, ), 0.0),
    ]
    if not output_all:
        decls += [
            L.ArrayDecl("double", J, (gdim * tdim, ), 0.0),
            L.ArrayDecl("double", detJ, (1, )),
            L.ArrayDecl("double", K, (tdim * gdim, ), 0.0),
        ]

    # By computing x and K at the cell midpoint once, we can optimize
    # the first iteration for each target point by initializing Xk = Xm
    # + Km * (xgoal - xm) which is the affine approximation starting at
    # the midpoint.
    midpoint_geometry = [
        L.Comment("Compute K = J^-1 and x at midpoint of cell"),
        L.ArrayDecl("double", xm, (gdim, ), 0.0),
        L.ArrayDecl("double", Km, (tdim * gdim, )),
        L.Call("compute_midpoint_geometry_{}".format(ir.name),
               (xm, J, coordinate_dofs)),
        L.Call("compute_jacobian_determinants_{}".format(ir.name),
               (detJ, one_point, J)),
        L.Call("compute_jacobian_inverses_{}".format(ir.name),
               (Km, one_point, J, detJ)),
    ]

    # declare xgoal = x[ip]
    # declare Xk = initial value
    newton_init = [
        L.Comment("Declare xgoal to hold the current x coordinate value"),
        L.ArrayDecl("const double", xgoal, (gdim, ), values=[x[ip][iv] for iv in range(gdim)]),
        L.Comment("Declare Xk iterate with initial value equal to reference cell midpoint"),
        L.ArrayDecl("double", Xk, (tdim, ), values=[Xm[c] for c in range(tdim)]),
        L.ForRanges(
            (j, 0, tdim), (i, 0, gdim),
            index_type=index_type,
            body=L.AssignAdd(Xk[j], Kmf[j, i] * (xgoal[i] - xm[i])))
    ]

    part1 = [
        L.Comment("Compute K = J^-1 for one point, (J and detJ are only used as"),
        L.Comment("intermediate storage inside compute_geometry, not used out here"),
        L.Call("compute_geometry_{}".format(ir.name),
               (xk, J, detJ, K, one_point, Xk, coordinate_dofs)),
    ]

    # Newton body with stopping criteria |dX|^2 < epsilon
    newton_body = part1 + [
        L.Comment("Declare dX increment to be computed, initialized to zero"),
        L.ArrayDecl("double", dX, (tdim, ), values=0.0),
        L.Comment("Compute dX[j] = sum_i K_ji * (x_i - x(Xk)_i)"),
        L.ForRanges(
            (j, 0, tdim), (i, 0, gdim),
            index_type=index_type,
            body=L.AssignAdd(dX[j], Kf[j, i] * (xgoal[i] - xk[i]))),
        L.Comment("Compute |dX|^2"),
        L.VariableDecl("double", dX2, value=0.0),
        L.ForRange(j, 0, tdim, index_type=index_type, body=L.AssignAdd(dX2, dX[j] * dX[j])),
        L.Comment("Break if converged (before X += dX such that X,J,detJ,K are consistent)"),
        L.If(L.LT(dX2, epsilon), L.Break()),
        L.Comment("Update Xk += dX"),
        L.ForRange(j, 0, tdim, index_type=index_type, body=L.AssignAdd(Xk[j], dX[j])),
    ]

    # Use this instead to have a fixed iteration number
    # alternative_newton_body = part1 + [
    #     L.Comment("Compute Xk[j] += sum_i K_ji * (x_i - x(Xk)_i)"),
    #     L.ForRanges(
    #         (j, 0, tdim), (i, 0, gdim),
    #         index_type=index_type,
    #         body=L.AssignAdd(Xk[j], Kf[j, i] * (xgoal[i] - xk[i]))),
    # ]

    # Carry out newton loop for each point
    point_loop = [
        L.ForRange(
            ip,
            0,
            num_points,
            index_type=index_type,
            body=[
                newton_init,
                # Loop until convergence
                L.ForRange(k, 0, max_iter, index_type=index_type, body=newton_body),
                # Copy X[ip] = Xk
                L.ForRange(j, 0, tdim, index_type=index_type, body=L.Assign(X[ip][j], Xk[j])),
            ])
    ]

    code = decls + midpoint_geometry + point_loop
    return code


def evaluate_reference_basis_derivatives_declaration(L, ir):
    scalar_coordinate_element_classname = ir.scalar_coordinate_finite_element_classname
    code = """
int evaluate_reference_basis_derivatives_{}(double* restrict reference_values,
    int order, int num_points, const double* restrict X);
""".format(scalar_coordinate_element_classname)
    return code


def compute_jacobians(L, ir):
    num_dofs = ir.num_scalar_coordinate_element_dofs
    scalar_coordinate_element_classname = ir.scalar_coordinate_finite_element_classname

    # Dimensions
    gdim = ir.geometric_dimension
    tdim = ir.topological_dimension
    num_points = L.Symbol("num_points")

    # Loop indices
    ip = L.Symbol("ip")
    i = L.Symbol("i")
    j = L.Symbol("j")
    d = L.Symbol("d")

    iz = L.Symbol("l")  # Array zeroing index

    # Input cell data
    coordinate_dofs = L.FlattenedArray(L.Symbol("coordinate_dofs"), dims=(num_dofs, gdim))

    # Output geometry
    J = L.FlattenedArray(L.Symbol("J"), dims=(num_points, gdim, tdim))

    # Input geometry
    X = L.FlattenedArray(L.Symbol("X"), dims=(num_points, tdim))

    # Computing table one point at a time instead of using num_points
    # will allow skipping dynamic allocation
    one_point = 1

    # Symbols for local basis derivatives table
    dphi_sym = L.Symbol("dphi")
    dphi = L.FlattenedArray(dphi_sym, dims=(num_dofs, tdim))

    # For each point, compute basis derivatives and accumulate into the right J
    code = [
        L.ArrayDecl("double", dphi_sym, (one_point * num_dofs * tdim, )),
        L.ForRange(
            iz, 0, num_points * gdim * tdim, index_type=index_type, body=L.Assign(J.array[iz],
                                                                                  0.0)),
        L.ForRange(
            ip,
            0,
            num_points,
            index_type=index_type,
            body=[
                L.Comment("Compute basis derivatives of coordinate element"),
                L.Call("evaluate_reference_basis_derivatives_{}".format(
                    scalar_coordinate_element_classname), (dphi_sym, 1, 1, L.AddressOf(X[ip, 0]))),
                L.Comment("Compute J"),
                L.ForRanges(
                    (i, 0, gdim), (j, 0, tdim), (d, 0, num_dofs),
                    index_type=index_type,
                    body=L.AssignAdd(J[ip, i, j], coordinate_dofs[d, i] * dphi[d, j]))
            ]),
    ]

    return code


def compute_jacobian_determinants(L, ir):
    # Dimensions
    gdim = ir.geometric_dimension
    tdim = ir.topological_dimension
    num_points = L.Symbol("num_points")

    # Loop indices
    ip = L.Symbol("ip")

    # Output geometry
    detJ = L.Symbol("detJ")[ip]

    # Input geometry
    J = L.FlattenedArray(L.Symbol("J"), dims=(num_points, gdim, tdim))

    # Assign det expression to detJ
    if gdim == tdim:
        body = L.Assign(detJ, det_nn(J[ip], gdim))
    elif tdim == 1:
        body = L.Assign(detJ, pdet_m1(L, J[ip], gdim))
    else:
        JTJ = L.Symbol("JTJ")
        body = [
            generate_compute_ATA(L, JTJ, J[ip], gdim, tdim),
            L.Assign(detJ, L.Sqrt(det_nn(JTJ, tdim))),
        ]

    # Carry out for all points
    code = L.ForRange(ip, 0, num_points, index_type=index_type, body=body)
    return code


def compute_jacobian_inverses(L, ir):
    # Dimensions
    gdim = ir.geometric_dimension
    tdim = ir.topological_dimension
    num_points = L.Symbol("num_points")

    # Loop indices
    ip = L.Symbol("ip")

    # Output geometry
    K = L.FlattenedArray(L.Symbol("K"), dims=(num_points, tdim, gdim))

    # Input geometry
    J = L.FlattenedArray(L.Symbol("J"), dims=(num_points, gdim, tdim))
    detJ = L.Symbol("detJ")

    # Assign to K[j][i] for each component j,i
    body = generate_assign_inverse(L, K[ip], J[ip], detJ[ip], gdim, tdim)

    # Carry out for all points
    return L.ForRange(ip, 0, num_points, index_type=index_type, body=body)


def compute_geometry(L, ir):
    # Class name
    classname = ir.name

    # Output geometry
    x = L.Symbol("x")
    J = L.Symbol("J")
    detJ = L.Symbol("detJ")
    K = L.Symbol("K")

    # Dimensions
    num_points = L.Symbol("num_points")

    # Input geometry
    X = L.Symbol("X")

    # Input cell data
    coordinate_dofs = L.Symbol("coordinate_dofs")

    # Just chain calls to other functions here
    code = [
        L.Call("compute_physical_coordinates_{}".format(classname),
               (x, num_points, X, coordinate_dofs)),
        L.Call("compute_jacobians_{}".format(classname), (J, num_points, X, coordinate_dofs)),
        L.Call("compute_jacobian_determinants_{}".format(classname),
               (detJ, num_points, J)),
        L.Call("compute_jacobian_inverses_{}".format(classname), (K, num_points, J, detJ)),
    ]

    return code


def compute_midpoint_geometry(L, ir):
    # Dimensions
    gdim = ir.geometric_dimension
    tdim = ir.topological_dimension
    num_dofs = ir.num_scalar_coordinate_element_dofs

    # Tables of coordinate basis function values and derivatives at
    # X=0 and X=midpoint available through ir. This is useful in
    # several geometry functions.
    tables = ir.tables

    # Check the table shapes against our expectations
    xm_table = tables["xm"]
    Jm_table = tables["Jm"]
    assert xm_table.shape == (num_dofs, )
    assert Jm_table.shape == (tdim, num_dofs)

    # TODO: Use epsilon parameter here?
    # TODO: Move to a more 'neutral' utility file
    from ffcx.ir.uflacs.elementtables import clamp_table_small_numbers
    xm_table = clamp_table_small_numbers(xm_table)
    Jm_table = clamp_table_small_numbers(Jm_table)

    # Table symbols
    phi_Xm = L.Symbol("phi_Xm")
    dphi_Xm = L.Symbol("dphi_Xm")

    # Table declarations
    table_decls = [
        L.ArrayDecl("const double", phi_Xm, sizes=xm_table.shape, values=xm_table),
        L.ArrayDecl("const double", dphi_Xm, sizes=Jm_table.shape, values=Jm_table),
    ]

    # Output geometry
    x = L.Symbol("x")
    J = L.Symbol("J")

    # Input cell data
    coordinate_dofs = L.Symbol("coordinate_dofs")
    coordinate_dofs = L.FlattenedArray(coordinate_dofs, dims=(num_dofs, gdim))
    Jf = L.FlattenedArray(J, dims=(gdim, tdim))

    i = L.Symbol("i")
    j = L.Symbol("j")
    d = L.Symbol("d")
    iz = L.Symbol("l")  # Index for zeroing arrays

    # Initialise arrays to zero
    init_array = [
        L.ForRange(iz, 0, gdim * tdim, index_type=index_type, body=L.Assign(Jf.array[iz], 0.0))
    ]

    xm_code = [
        L.Comment("Compute x"),
        L.ForRange(
            i,
            0,
            gdim,
            index_type=index_type,
            body=[
                L.Assign(x[i], 0.0),
                L.ForRange(
                    d,
                    0,
                    num_dofs,
                    index_type=index_type,
                    body=L.AssignAdd(x[i], coordinate_dofs[d, i] * phi_Xm[d]))
            ]),
    ]

    Jm_code = [
        L.Comment("Compute J"),
        L.ForRanges(
            (i, 0, gdim), (j, 0, tdim), (d, 0, num_dofs),
            index_type=index_type,
            body=L.AssignAdd(Jf[i, j], coordinate_dofs[d, i] * dphi_Xm[j, d])),
    ]

    # Reuse functions for detJ and K
    code = table_decls + init_array + xm_code + Jm_code
    return code


def generator(ir, parameters):
    """Generate UFC code for a coordinate mapping."""

    d = {}

    # Attributes
    d["factory_name"] = ir.name
    d["signature"] = "\"{}\"".format(ir.signature)
    d["geometric_dimension"] = ir.geometric_dimension
    d["topological_dimension"] = ir.topological_dimension
    d["cell_shape"] = ir.cell_shape

    import ffcx.codegeneration.C.cnodes as L

    statements = compute_physical_coordinates(L, ir)
    d["compute_physical_coordinates"] = L.StatementList(statements)
    d["evaluate_reference_basis_declaration"] = evaluate_reference_basis_declaration(L, ir)

    statements = compute_reference_coordinates(L, ir)
    d["compute_reference_coordinates"] = L.StatementList(statements)

    statements = compute_reference_geometry(L, ir)
    d["compute_reference_geometry"] = L.StatementList(statements)
    d["evaluate_reference_basis_derivatives_declaration"] = evaluate_reference_basis_derivatives_declaration(L, ir)

    statements = compute_jacobians(L, ir)
    d["compute_jacobians"] = L.StatementList(statements)

    d["compute_jacobian_determinants"] = compute_jacobian_determinants(L, ir)
    d["compute_jacobian_inverses"] = compute_jacobian_inverses(L, ir)

    statements = compute_geometry(L, ir)
    d["compute_geometry"] = L.StatementList(statements)

    statements = compute_midpoint_geometry(L, ir)
    d["compute_midpoint_geometry"] = L.StatementList(statements)

    # Check that no keys are redundant or have been missed
    from string import Formatter
    fields = [
        fname for _, fname, _, _ in Formatter().parse(ufc_coordinate_mapping.factory) if fname
    ]
    assert set(fields) == set(
        d.keys()), "Mismatch between keys in template and in formattting dict."

    # Format implementation code
    implementation = ufc_coordinate_mapping.factory.format_map(d)

    # Format declaration
    declaration = ufc_coordinate_mapping.declaration.format(factory_name=ir.name)

    return declaration, implementation
