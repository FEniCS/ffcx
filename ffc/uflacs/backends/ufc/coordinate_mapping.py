# -*- coding: utf-8 -*-
# Copyright (C) 2015-2017 Martin Sandve Alnæs
#
# This file is part of UFLACS.
#
# UFLACS is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# UFLACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with UFLACS. If not, see <http://www.gnu.org/licenses/>.

from ffc.uflacs.backends.ufc.generator import ufc_generator
from ffc.uflacs.backends.ufc.utils import generate_return_new


# TODO: Test everything here! Cover all combinations of gdim,tdim=1,2,3!


index_type = "std::size_t"


### Code generation utilities:

def generate_compute_ATA(L, ATA, A, m, n, index_prefix=""):
    "Generate code to declare and compute ATA[i,j] = sum_k A[k,i]*A[k,j] with given A shaped (m,n)."
    # Loop indices
    i = L.Symbol(index_prefix + "i")
    j = L.Symbol(index_prefix + "j")
    k = L.Symbol(index_prefix + "k")

    # Build A^T*A matrix
    code = [
        L.ArrayDecl("double", ATA, (n, n), values=0),
        L.ForRanges(
            (i, 0, n),
            (j, 0, n),
            (k, 0, m),
            index_type=index_type,
            body=L.AssignAdd(ATA[i, j], A[k, i] * A[k, j])
        ),
    ]
    return L.StatementList(code)

def cross_expr(a, b):
    def cr(i, j):
        return a[i]*b[j] - a[j]*b[i]
    return [cr(1, 2), cr(2, 0), cr(0, 1)]

def generate_cross_decl(c, a, b):
    return L.ArrayDecl("double", c, values=cross_expr(a, b))

### Inline math expressions:

def det_22(B, i, j, k, l):
    return B[i, k]*B[j, l] - B[i, l]*B[j, k]

def codet_nn(A, rows, cols):
    n = len(rows)
    if n == 2:
        return det_22(A, rows[0], rows[1], cols[0], cols[1])
    else:
        r = rows[0]
        subrows = rows[1:]
        parts = []
        for i, c in enumerate(cols):
            subcols = cols[i+1:] + cols[:i]
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
    A2 = A[0,0]*A[0,0]
    for i in range(1, m):
        A2 = A2 + A[i,0]*A[i,0]
    return L.Sqrt(A2)

def adj_expr_2x2(A):
    return [[ A[1, 1], -A[0, 1]],
            [-A[1, 0],  A[0, 0]]]

def adj_expr_3x3(A):
    return [[A[2, 2]*A[1, 1] - A[1, 2]*A[2, 1],
             A[0, 2]*A[2, 1] - A[0, 1]*A[2, 2],
             A[0, 1]*A[1, 2] - A[0, 2]*A[1, 1]],
            [A[1, 2]*A[2, 0] - A[2, 2]*A[1, 0],
             A[2, 2]*A[0, 0] - A[0, 2]*A[2, 0],
             A[0, 2]*A[1, 0] - A[1, 2]*A[0, 0]],
            [A[1, 0]*A[2, 1] - A[2, 0]*A[1, 1],
             A[0, 1]*A[2, 0] - A[0, 0]*A[2, 1],
             A[0, 0]*A[1, 1] - A[0, 1]*A[1, 0]]]

def generate_assign_inverse(L, K, J, detJ, gdim, tdim):
    if gdim == tdim:
        if gdim == 1:
            return L.Assign(K[0, 0], L.LiteralFloat(1.0) / J[0, 0])
        elif gdim == 2:
            adj_values = adj_expr_2x2(J)
        elif gdim == 3:
            adj_values = adj_expr_3x3(J)
        else:
            error("Not implemented.")
        return L.StatementList([L.Assign(K[j,i], adj_values[j][i] / detJ)
                                for j in range(tdim) for i in range(gdim)])
    else:
        if tdim == 1:
            # Simpler special case for embedded 1d
            prods = [J[i,0]*J[i,0] for i in range(gdim)]
            s = sum(prods[1:], prods[0])
            return L.StatementList([L.Assign(K[0,i], J[i,0] / s)
                                    for i in range(gdim)])
        else:
            # Generic formulation of Penrose-Moore pseudo-inverse of J: (J.T*J)^-1 * J.T
            i = L.Symbol("i")
            j = L.Symbol("j")
            k = L.Symbol("k")
            JTJ = L.Symbol("JTJ")
            JTJf = L.FlattenedArray(JTJ, dims=(tdim, tdim))
            detJTJ = detJ*detJ
            JTJinv = L.Symbol("JTJinv")
            JTJinvf = L.FlattenedArray(JTJinv, dims=(tdim, tdim))
            code = [
                L.Comment("Compute J^T J"),
                L.ArrayDecl("double", JTJ, (tdim*tdim,), values=0),
                L.ForRanges(
                    (k, 0, tdim),
                    (j, 0, tdim),
                    (i, 0, gdim),
                    index_type=index_type,
                    body=L.AssignAdd(JTJf[k, j], J[i, k] * J[i, j])
                ),
                L.Comment("Compute inverse(J^T J)"),
                L.ArrayDecl("double", JTJinv, (tdim*tdim,)),
                generate_assign_inverse(L, JTJinvf, JTJf, detJTJ, tdim, tdim),
                L.Comment("Compute K = inverse(J^T J) * J"),
                L.ForRanges(
                    (k, 0, tdim),
                    (i, 0, gdim),
                    index_type=index_type,
                    body=L.Assign(K[k, i], L.LiteralFloat(0.0))
                ),
                L.ForRanges(
                    (k, 0, tdim),
                    (i, 0, gdim),
                    (j, 0, tdim),
                    index_type=index_type,
                    body=L.AssignAdd(K[k, i], JTJinvf[k, j] * J[i, j])
                ),
                ]
            return L.StatementList(code)


class ufc_coordinate_mapping(ufc_generator):
    "Each function maps to a keyword in the template. See documentation of ufc_generator."
    def __init__(self):
        ufc_generator.__init__(self, "coordinate_mapping")

    def cell_shape(self, L, ir):
        name = ir["cell_shape"]
        return L.Return(L.Symbol("ufc::shape::" + name))

    def topological_dimension(self, L, topological_dimension):
        return L.Return(topological_dimension)

    def geometric_dimension(self, L, geometric_dimension):
        return L.Return(geometric_dimension)

    def create_coordinate_finite_element(self, L, ir):
        classname = ir["create_coordinate_finite_element"]
        return generate_return_new(L, classname, factory=ir["jit"])

    def create_coordinate_dofmap(self, L, ir):
        classname = ir["create_coordinate_dofmap"]
        return generate_return_new(L, classname, factory=ir["jit"])

    def compute_physical_coordinates(self, L, ir):
        num_dofs = ir["num_scalar_coordinate_element_dofs"]
        scalar_coordinate_element_classname = ir["scalar_coordinate_finite_element_classname"]

        # Dimensions
        gdim = ir["geometric_dimension"]
        tdim = ir["topological_dimension"]
        num_points = L.Symbol("num_points")

        # Loop indices
        ip = L.Symbol("ip")
        i = L.Symbol("i")
        j = L.Symbol("j")
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
        phi = L.FlattenedArray(phi_sym, dims=(num_dofs,))

        # For each point, compute basis values and accumulate into the right x
        code = [
            # Define scalar finite element instance
            # (stateless, so placing this on the stack is basically free)
            L.VariableDecl(scalar_coordinate_element_classname, "xelement"),
            L.ArrayDecl("double", phi_sym, (one_point*num_dofs,)),
            L.ForRange(ip, 0, num_points, index_type=index_type, body=[
                L.Comment("Compute basis values of coordinate element"),
                L.Call("xelement.evaluate_reference_basis", (phi_sym, 1, L.AddressOf(X[ip, 0]))),
                L.Comment("Compute x"),
                L.ForRanges(
                    (i, 0, gdim),
                    (d, 0, num_dofs),
                    index_type=index_type,
                    body=L.AssignAdd(x[ip, i], coordinate_dofs[d, i]*phi[d])
                )
            ])
        ]
        return code

    def compute_reference_geometry(self, L, ir):
        degree = ir["coordinate_element_degree"]
        if degree == 1:
            # Special case optimized for affine mesh (possibly room for further optimization)
            return self._compute_reference_coordinates_affine(L, ir, output_all=True)
        else:
            # General case with newton loop to solve F(X) = x(X) - x0 = 0
            return self._compute_reference_coordinates_newton(L, ir, output_all=True)

    # TODO: Maybe we don't need this version, see what we need in dolfin first
    def compute_reference_coordinates(self, L, ir):
        degree = ir["coordinate_element_degree"]
        if degree == 1:
            # Special case optimized for affine mesh (possibly room for further optimization)
            return self._compute_reference_coordinates_affine(L, ir)
        else:
            # General case with newton loop to solve F(X) = x(X) - x0 = 0
            return self._compute_reference_coordinates_newton(L, ir)

    def _compute_reference_coordinates_affine(self, L, ir, output_all=False):
        # Dimensions
        gdim = ir["geometric_dimension"]
        tdim = ir["topological_dimension"]
        cellname = ir["cell_shape"]
        num_points = L.Symbol("num_points")

        # Number of dofs for a scalar component
        num_dofs = ir["num_scalar_coordinate_element_dofs"]

        # Loop indices
        ip = L.Symbol("ip") # point
        i = L.Symbol("i")   # gdim
        j = L.Symbol("j")   # tdim
        k = L.Symbol("k")   # sum iteration

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
        cell_orientation = L.Symbol("cell_orientation")

        if output_all:
            decls = []
        else:
            decls = [
                L.ArrayDecl("double", Jsym, sizes=(gdim*tdim,)),
                L.ArrayDecl("double", detJsym, sizes=(1,)),
                L.ArrayDecl("double", Ksym, sizes=(tdim*gdim,)),
            ]

        # Tables of coordinate basis function values and derivatives at
        # X=0 and X=midpoint available through ir. This is useful in
        # several geometry functions.
        tables = ir["tables"]

        # Check the table shapes against our expectations
        x_table = tables["x0"]
        J_table = tables["J0"]
        assert x_table.shape == (num_dofs,)
        assert J_table.shape == (tdim, num_dofs)

        # TODO: Use epsilon parameter here?
        # TODO: Move to a more 'neutral' utility file
        from ffc.uflacs.elementtables import clamp_table_small_numbers
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
            L.ArrayDecl("double", x0, sizes=(gdim,), values=0),
            L.ForRanges(
                (i, 0, gdim),
                (k, 0, num_dofs),
                index_type=index_type,
                body=L.AssignAdd(x0[i], coordinate_dofs[k, i] * phi_X0[k])
            ),
        ]

        # For more convenient indexing
        detJ = detJsym[0]
        J = L.FlattenedArray(Jsym, dims=(gdim, tdim))
        K = L.FlattenedArray(Ksym, dims=(tdim, gdim))

        # Compute J = J(X=0) (optimized by precomputing basis at X=0)
        compute_J0 = [
            L.ForRanges(
                (i, 0, gdim),
                (j, 0, tdim),
                index_type=index_type,
                body=[
                    L.Assign(J[i, j], 0.0),
                    L.ForRange(k, 0, num_dofs, index_type=index_type, body=
                        L.AssignAdd(J[i, j], coordinate_dofs[k, i] * dphi_X0[j, k])
                    )
                ]
            ),
            ]

        # Compute K = inv(J) (and intermediate value det(J))
        compute_K0 = [
            L.Call("compute_jacobian_determinants", (detJsym, 1, Jsym, cell_orientation)),
            L.Call("compute_jacobian_inverses", (Ksym, 1, Jsym, detJsym)),
            ]

        # Compute X = K0*(x-x0) for each physical point x
        compute_X = [
            L.ForRanges(
                (ip, 0, num_points),
                (j, 0, tdim),
                (i, 0, gdim),
                index_type=index_type,
                body=L.AssignAdd(X[ip, j], K[j, i]*(x[ip, i] - x0[i]))
            )
            ]

        # Stitch it together
        code = table_decls + decls + compute_x0 + compute_J0 + compute_K0 + compute_X
        return code

    def _compute_reference_coordinates_newton(self, L, ir, output_all=False):
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
        gdim = ir["geometric_dimension"]
        tdim = ir["topological_dimension"]
        cellname = ir["cell_shape"]
        num_points = L.Symbol("num_points")

        degree = ir["coordinate_element_degree"]

        # Computing table one point at a time instead of vectorized
        # over num_points will allow skipping dynamic allocation
        one_point = 1

        # Loop indices
        ip = L.Symbol("ip") # point
        i = L.Symbol("i")   # gdim
        j = L.Symbol("j")   # tdim
        k = L.Symbol("k")   # iteration

        # Input cell data
        coordinate_dofs = L.Symbol("coordinate_dofs")
        cell_orientation = L.Symbol("cell_orientation")

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
        #max_iter = L.Symbol("iterations")
        #epsilon = L.Symbol("epsilon")
        dX2 = L.Symbol("dX2")

        # Wrap K as flattened array for convenient indexing Kf[j,i]
        Kf = L.FlattenedArray(K, dims=(tdim, gdim))

        decls = [
            L.Comment("Declare intermediate arrays to hold results of compute_geometry call"),
            L.ArrayDecl("double", xk, (gdim,), 0.0),
            ]
        if not output_all:
            decls += [
                L.ArrayDecl("double", J, (gdim*tdim,), 0.0),
                L.ArrayDecl("double", detJ, (1,)),
                L.ArrayDecl("double", K, (tdim*gdim,), 0.0),
                ]

        # By computing x and K at the cell midpoint once,
        # we can optimize the first iteration for each target point
        # by initializing Xk = Xm + Km * (xgoal - xm)
        # which is the affine approximation starting at the midpoint.
        midpoint_geometry = [
            L.Comment("Compute K = J^-1 and x at midpoint of cell"),
            L.ArrayDecl("double", xm, (gdim,), 0.0),
            L.ArrayDecl("double", Km, (tdim*gdim,)),
            L.Call("compute_midpoint_geometry", (xm, J, coordinate_dofs)),
            L.Call("compute_jacobian_determinants", (detJ, one_point, J, cell_orientation)),
            L.Call("compute_jacobian_inverses", (Km, one_point, J, detJ)),
        ]

        # declare xgoal = x[ip]
        # declare Xk = initial value
        newton_init = [
            L.Comment("Declare xgoal to hold the current x coordinate value"),
            L.ArrayDecl("const double", xgoal, (gdim,), values=[x[ip][iv] for iv in range(gdim)]),
            L.Comment("Declare Xk iterate with initial value equal to reference cell midpoint"),
            L.ArrayDecl("double", Xk, (tdim,), values=[Xm[c] for c in range(tdim)]),
            L.ForRanges(
                (j, 0, tdim),
                (i, 0, gdim),
                index_type=index_type,
                body=L.AssignAdd(Xk[j], Kf[j,i] * (xgoal[i] - xm[i]))
            )
        ]

        part1 = [
            L.Comment("Compute K = J^-1 for one point, (J and detJ are only used as"),
            L.Comment("intermediate storage inside compute_geometry, not used out here"),
            L.Call("compute_geometry",
                   (xk, J, detJ, K, one_point,
                    Xk, coordinate_dofs, cell_orientation)),
        ]

        # Newton body with stopping criteria |dX|^2 < epsilon
        newton_body = part1 + [
            L.Comment("Declare dX increment to be computed, initialized to zero"),
            L.ArrayDecl("double", dX, (tdim,), values=0.0),

            L.Comment("Compute dX[j] = sum_i K_ji * (x_i - x(Xk)_i)"),
            L.ForRanges(
                (j, 0, tdim),
                (i, 0, gdim),
                index_type=index_type,
                body=L.AssignAdd(dX[j], Kf[j,i]*(xgoal[i] - xk[i]))
            ),

            L.Comment("Compute |dX|^2"),
            L.VariableDecl("double", dX2, value=0.0),
            L.ForRange(j, 0, tdim, index_type=index_type,
                body=L.AssignAdd(dX2, dX[j]*dX[j])),

            L.Comment("Break if converged (before X += dX such that X,J,detJ,K are consistent)"),
            L.If(L.LT(dX2, epsilon), L.Break()),

            L.Comment("Update Xk += dX"),
            L.ForRange(j, 0, tdim, index_type=index_type,
                body=L.AssignAdd(Xk[j], dX[j])),
        ]

        # Use this instead to have a fixed iteration number
        alternative_newton_body = part1 + [
            L.Comment("Compute Xk[j] += sum_i K_ji * (x_i - x(Xk)_i)"),
            L.ForRanges(
                (j, 0, tdim),
                (i, 0, gdim),
                index_type=index_type,
                body=L.AssignAdd(Xk[j], Kf[j,i]*(xgoal[i] - xk[i]))
            ),
        ]

        # Carry out newton loop for each point
        point_loop = [
            L.ForRange(ip, 0, num_points, index_type=index_type, body=[
                newton_init,
                # Loop until convergence
                L.ForRange(k, 0, max_iter, index_type=index_type,
                    body=newton_body),
                # Copy X[ip] = Xk
                L.ForRange(j, 0, tdim, index_type=index_type,
                    body=L.Assign(X[ip][j], Xk[j])),
            ])
        ]

        code = decls + midpoint_geometry + point_loop
        return code

    def compute_jacobians(self, L, ir):
        num_dofs = ir["num_scalar_coordinate_element_dofs"]
        scalar_coordinate_element_classname = ir["scalar_coordinate_finite_element_classname"]

        # Dimensions
        gdim = ir["geometric_dimension"]
        tdim = ir["topological_dimension"]
        num_points = L.Symbol("num_points")

        # Loop indices
        ip = L.Symbol("ip")
        i = L.Symbol("i")
        j = L.Symbol("j")
        d = L.Symbol("d")

        # Input cell data
        coordinate_dofs = L.FlattenedArray(L.Symbol("coordinate_dofs"), dims=(num_dofs, gdim))

        # Output geometry
        J = L.FlattenedArray(L.Symbol("J"), dims=(num_points, gdim, tdim))

        # Input geometry
        X = L.FlattenedArray(L.Symbol("X"), dims=(num_points, tdim))

        # Computing table one point at a time instead of using
        # num_points will allow skipping dynamic allocation
        one_point = 1

        # Symbols for local basis derivatives table
        dphi_sym = L.Symbol("dphi")
        dphi = L.FlattenedArray(dphi_sym, dims=(num_dofs, tdim))

        # For each point, compute basis derivatives and accumulate into the right J
        code = [
            L.VariableDecl(scalar_coordinate_element_classname, "xelement"),
            L.ArrayDecl("double", dphi_sym, (one_point*num_dofs*tdim,)),
            L.ForRange(ip, 0, num_points, index_type=index_type, body=[
                L.Comment("Compute basis derivatives of coordinate element"),
                L.Call("xelement.evaluate_reference_basis_derivatives",
                       (dphi_sym, 1, 1, L.AddressOf(X[ip, 0]))),
                L.Comment("Compute J"),
                L.ForRanges(
                    (i, 0, gdim),
                    (j, 0, tdim),
                    (d, 0, num_dofs),
                    index_type=index_type,
                    body=L.AssignAdd(J[ip, i, j], coordinate_dofs[d, i]*dphi[d, j])
                )
            ])
        ]
        return code

    def compute_jacobian_determinants(self, L, ir):
        # Dimensions
        gdim = ir["geometric_dimension"]
        tdim = ir["topological_dimension"]
        num_points = L.Symbol("num_points")

        # Loop indices
        ip = L.Symbol("ip")

        # Output geometry
        detJ = L.Symbol("detJ")[ip]

        # Input geometry
        J = L.FlattenedArray(L.Symbol("J"), dims=(num_points, gdim, tdim))
        cell_orientation = L.Symbol("cell_orientation")
        orientation_scaling = L.Conditional(L.EQ(cell_orientation, 1), -1.0, +1.0)

        # Assign det expression to detJ
        if gdim == tdim:
            body = L.Assign(detJ, det_nn(J[ip], gdim))
        elif tdim == 1:
            body = L.Assign(detJ, orientation_scaling*pdet_m1(L, J[ip], gdim))
        #elif tdim == 2 and gdim == 3:
        #    body = L.Assign(detJ, orientation_scaling*pdet_32(L, J[ip])) # Possible optimization not implemented here
        else:
            JTJ = L.Symbol("JTJ")
            body = [
                generate_compute_ATA(L, JTJ, J[ip], gdim, tdim),
                L.Assign(detJ, orientation_scaling*L.Sqrt(det_nn(JTJ, tdim))),
                ]

        # Carry out for all points
        code = L.ForRange(ip, 0, num_points, index_type=index_type, body=body)
        return code

    def compute_jacobian_inverses(self, L, ir):
        # Dimensions
        gdim = ir["geometric_dimension"]
        tdim = ir["topological_dimension"]
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

    def compute_geometry(self, L, ir):
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
        cell_orientation = L.Symbol("cell_orientation")

        # All arguments
        args = (x, J, detJ, K, num_points, X, coordinate_dofs, cell_orientation)

        # Just chain calls to other functions here
        code = [
            L.Call("compute_physical_coordinates", (x, num_points, X, coordinate_dofs)),
            L.Call("compute_jacobians", (J, num_points, X, coordinate_dofs)),
            L.Call("compute_jacobian_determinants", (detJ, num_points, J, cell_orientation)),
            L.Call("compute_jacobian_inverses", (K, num_points, J, detJ)),
            ]
        return code

    def compute_midpoint_geometry(self, L, ir):
        # Dimensions
        gdim = ir["geometric_dimension"]
        tdim = ir["topological_dimension"]
        num_dofs = ir["num_scalar_coordinate_element_dofs"]

        # Tables of coordinate basis function values and derivatives at
        # X=0 and X=midpoint available through ir. This is useful in
        # several geometry functions.
        tables = ir["tables"]

        # Check the table shapes against our expectations
        xm_table = tables["xm"]
        Jm_table = tables["Jm"]
        assert xm_table.shape == (num_dofs,)
        assert Jm_table.shape == (tdim, num_dofs)

        # TODO: Use epsilon parameter here?
        # TODO: Move to a more 'neutral' utility file
        from ffc.uflacs.elementtables import clamp_table_small_numbers
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

        # Symbol for ufc_geometry cell midpoint definition
        cellname = ir["cell_shape"]
        Xm = L.Symbol("%s_midpoint" % cellname)

        # Output geometry
        x = L.Symbol("x")
        J = L.Symbol("J")
        detJ = L.Symbol("detJ")
        K = L.Symbol("K")

        # Dimensions
        num_points = 1

        # Input geometry
        X = Xm

        # Input cell data
        coordinate_dofs = L.Symbol("coordinate_dofs")
        coordinate_dofs = L.FlattenedArray(coordinate_dofs, dims=(num_dofs, gdim))
        Jf = L.FlattenedArray(J, dims=(gdim, tdim))

        i = L.Symbol("i")
        j = L.Symbol("j")
        d = L.Symbol("d")

        xm_code = [
            L.Comment("Compute x"),
            L.ForRanges(
                (i, 0, gdim),
                (d, 0, num_dofs),
                index_type=index_type,
                body=L.AssignAdd(x[i], coordinate_dofs[d, i]*phi_Xm[d])
            ),
        ]

        Jm_code = [
            L.Comment("Compute J"),
            L.ForRanges(
                (i, 0, gdim),
                (j, 0, tdim),
                (d, 0, num_dofs),
                index_type=index_type,
                body=L.AssignAdd(Jf[i, j], coordinate_dofs[d, i]*dphi_Xm[j, d])
            ),
        ]

        # Reuse functions for detJ and K
        code = table_decls + xm_code + Jm_code
        return code
