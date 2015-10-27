# -*- coding: utf-8 -*-
# Copyright (C) 2015-2015 Martin Sandve Alnæs
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

from uflacs.backends.ufc.generator import ufc_generator

### Code generation utilities:

def generate_compute_ATA(L, ATA, A, m, n, index_prefix=""):
    "Generate code to declare and compute ATA[i,j] = sum_k A[k,i]*A[k,j] with given A shaped (m,n)."
    # Loop indices
    i = L.Symbol(index_prefix + "i")
    j = L.Symbol(index_prefix + "j")
    k = L.Symbol(index_prefix + "k")

    # Build A^T*A matrix
    code = [
        L.ArrayDecl("double", ATA, sizes=(n, n), values=0),
        L.ForRange(i, 0, n, body=
            L.ForRange(j, 0, n, body=
                L.ForRange(k, 0, m, body=
                    L.AssignAdd(ATA[i, j], A[k, i] * A[k, j])))),
        ]
    return L.StatementList(code)

def generate_accumulation_loop(dst, expr, indices):
    body = L.AssignAdd(dst, expr)
    for index, begin, end in reversed(indices):
        body = L.ForRange(index, begin, end, body=body)
    return body

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
    return L.Call("sqrt", A2)

class ufc_coordinate_mapping(ufc_generator):
    def __init__(self):
        ufc_generator.__init__(self, "coordinate_mapping")

    def cell_shape(self, L, ir):
        name = ir["cell_shape"]
        return L.Return(L.Symbol(name))

    def topological_dimension(self, L, ir):
        "Default implementation of returning topological dimension fetched from ir."
        value = ir["topological_dimension"]
        return L.Return(L.LiteralInt(value))

    def geometric_dimension(self, L, ir):
        "Default implementation of returning geometric dimension fetched from ir."
        value = ir["geometric_dimension"]
        return L.Return(L.LiteralInt(value))

    def create_coordinate_finite_element(self, L, ir):
        classname = ir["create_coordinate_finite_element"] # FIXME: ffc passes class id not name
        return L.Return(L.New(classname))

    def create_coordinate_dofmap(self, L, ir):
        classname = ir["create_coordinate_dofmap"] # FIXME: ffc passes class id not name
        return L.Return(L.New(classname))

    def compute_physical_coordinates(self, L, ir): # FIXME: Fix jacobian implementation first then mirror solutions here.
        if 1: return L.Comment("FIXME")

        # Dimensions
        gdim = ir["geometric_dimension"]
        tdim = ir["topological_dimension"]
        num_points = L.Symbol("num_points")

        # Loop indices
        ip = L.Symbol("ip")
        i = L.Symbol("i")
        j = L.Symbol("j")
        dof = L.Symbol("dof")

        # Input cell data
        coordinate_dofs = L.Symbol("coordinate_dofs")
        cell_orientation = L.Symbol("cell_orientation")

        # Output geometry
        x = L.FlattenedArray(L.Symbol("x"), dims=(num_points, gdim))

        # Input geometry
        X = L.FlattenedArray(L.Symbol("X"), dims=(num_points, tdim))

        # FIXME: Almost exactly like jacobian implementation, solve issues there first
        num_dofs = 3
        phi = L.Symbol("phi") # FIXME: Compute basis functions in reference coordinates

        body = L.ForRange(i, 0, gdim, body=L.StatementList([
            # Assign to x[ip][i]
            L.Assign(x[ip][i], 0),
            L.ForRange(dof, 0, num_dofs, body=
                L.AssignAdd(x[i], coordinate_dofs[dof*gdim + i] * phi[ip,dof])
                )
            ]))

        # Carry out for all points
        return L.ForRange(ip, 0, num_points, body=body)

    def compute_reference_coordinates(self, L, ir):
        degree = ir["coordinate_element_degree"]
        if degree == 1:
            # Special case optimized for affine mesh (possibly room for further optimization)
            return self._compute_reference_coordinates_affine(L, ir)
        else:
            # General case with newton loop to solve F(X) = x(X) - x0 = 0
            return self._compute_reference_coordinates_newton(L, ir)

    def _compute_reference_coordinates_affine(self, L, ir):
        # Dimensions
        gdim = ir["geometric_dimension"]
        tdim = ir["topological_dimension"]
        cellname = ir["cell_shape"]
        num_points = L.Symbol("num_points")

        # Loop indices
        ip = L.Symbol("ip") # point
        i = L.Symbol("i")   # gdim
        j = L.Symbol("j")   # tdim
        k = L.Symbol("k")   # sum iteration

        # Output geometry
        X = L.FlattenedArray(L.Symbol("X"), dims=(num_points, tdim))

        # Input geometry
        x = L.FlattenedArray(L.Symbol("x"), dims=(num_points, gdim))

        # Input cell data
        coordinate_dofs = L.Symbol("coordinate_dofs")

        # Tables of coordinate basis function values and derivatives at
        # X=0 and X=midpoint available through ir. This is useful in
        # several geometry functions.
        tables = ir["tables"]

        # Table symbols
        phi_X0 = L.Symbol("phi_X0")
        dphi_X0 = L.Symbol("dphi_X0")

        # Interpret and check the table shapes
        # FIXME: Massage tables in ffc such that 1 point-dim is gone
        x_table = tables["x0"]
        J_table = tables["J0"]
        assert len(x_table.shape)
        num_dofs, = x_table.shape # Number of dofs for a scalar component FIXME: Add to ir and compare here instead
        assert J_table.shape == (tdim,) + x_table.shape

        # Table declarations
        table_decls = [
            L.ArrayDecl("static const double", phi_X0, sizes=x_table.shape, values=tables["x0"]),
            L.ArrayDecl("static const double", dphi_X0, sizes=J_table.shape, values=tables["J0"]),
            ]

        # This is the assumed shape of coordinate_dofs below
        xdofs = L.FlattenedArray(coordinate_dofs, dims=(num_dofs, gdim))

        # Compute x0 = x(X=0) (optimized by precomputing basis at X=0)
        x0 = L.Symbol("x0")
        compute_x0 = [
            L.ArrayDecl("double", x0, sizes=(gdim,), values=0),
            L.ForRange(i, 0, gdim, body=
                L.ForRange(k, 0, num_dofs, body=
                    L.AssignAdd(x0[i], xdofs[k, i] * phi_X0[k]))),
            ]

        # Compute J0 = J(X=0) (optimized by precomputing basis at X=0)
        J0 = L.Symbol("J0")
        compute_J0 = [
            L.ArrayDecl("double", J0, sizes=(gdim*tdim,), values=0),
            L.ForRange(i, 0, gdim, body=
                L.ForRange(j, 0, tdim, body=
                    L.ForRange(k, 0, num_dofs, body=
                        L.AssignAdd(J0[i*tdim + j], xdofs[k, i] * dphi_X0[j, k])))),
            ]

        # Compute K0 = inv(J0) (and intermediate value det(J0))
        detJ0 = L.Symbol("detJ0")
        K0 = L.Symbol("K0")
        compute_K0 = [
            L.ArrayDecl("double", detJ0, sizes=(1,)),
            L.Call("compute_jacobian_determinants", (detJ0, 1, J0)),
            L.ArrayDecl("double", K0, sizes=(tdim*gdim,)),
            L.Call("compute_jacobian_inverses", (K0, 1, J0, detJ0)),
            ]

        # Compute X = K0*(x-x0) for each physical point x
        compute_X = [
            L.ForRange(ip, 0, num_points, body=
                L.ForRange(j, 0, tdim, body=
                    L.ForRange(i, 0, gdim, body=
                        L.AssignAdd(X[ip, j], K0[j, i]*(x[ip, i] - x0[i]))))),
            ]

        # Stitch it together
        code = L.StatementList(table_decls + compute_x0 + compute_J0 + compute_K0 + compute_X)
        return code

    def _compute_reference_coordinates_newton(self, L, ir): # FIXME: Get degree of mapping. Otherwise looks mostly ok on inspection. TODO: Test and determine stopping criteria to use.
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

        # TODO: Implement:
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

        # Symbol for ufc_geometry cell midpoint definition
        mp = L.Symbol("%s_midpoint" % cellname)

        dX2 = L.Symbol("dX2")
        epsilon = L.LiteralFloat(1e-14) # L.Symbol("epsilon") # TODO: Choose good convergence criteria

        # Wrap K as flattened array for convenient indexing Kf[j,i]
        Kf = L.FlattenedArray(K, dims=(tdim, gdim))

        # declare xgoal = x[ip]
        # declare Xk = initial value
        newton_init = [
            L.Comment("Declare xgoal to hold the current x coordinate value"),
            L.ArrayDecl("const double", xgoal, (gdim,), values=[x[ip][iv] for iv in range(gdim)]),
            L.Comment("Declare Xk iterate with initial value equal to reference cell midpoint"),
            L.ArrayDecl("double", Xk, (tdim,), values=[mp[c] for c in range(tdim)]),
            ]

        part1 = [
            L.Comment("Declare intermediate arrays to hold results of compute_geometry call"),
            L.ArrayDecl("double", xk, (gdim,)),
            L.ArrayDecl("double", J, (gdim*tdim,)),
            L.ArrayDecl("double", detJ, (1,)),
            L.ArrayDecl("double", K, (tdim*gdim,)),

            L.Comment("Compute K = J^-1 for one point, (J and detJ are only used as"),
            L.Comment("intermediate storage inside compute_geometry, not used out here"),
            L.Call("compute_geometry",
                   (xk, J, detJ, K, one_point,
                    Xk, coordinate_dofs, cell_orientation)),
            ]

        part2a = [
            L.Comment("Declare dX increment to be computed, initialized to zero"),
            L.ArrayDecl("double", dX, (tdim,), values=0.0),

            L.Comment("Compute dX[j] = sum_i K_ji * (x(Xk)_i - x_i)"),
            L.ForRange(j, 0, tdim, body=
                       L.ForRange(i, 0, gdim, body=
                                  L.AssignAdd(dX[j], Kf[j,i]*(xk[i] - xgoal[i])))),

            L.Comment("Update Xk -= dX"),
            L.ForRange(j, 0, tdim, body=L.AssignSub(Xk[j], dX[j])),

            # TODO: If we set epsilon strict and break _before_ Xk -= dX, we can output consistent (J, detJ, K) together with X
            L.Comment("Compute |dX|^2"),
            L.VariableDecl("double", dX2, value=0.0),
            L.ForRange(j, 0, tdim, body=L.AssignAdd(dX2, dX[j]*dX[j])),

            L.Comment("Break if converged"),
            L.If(L.LT(dX2, epsilon), L.Break()),
            ]

        part2b = [
            L.Comment("Compute Xk[j] -= sum_i K_ji * (x(Xk)_i - x_i)"),
            L.ForRange(j, 0, tdim, body=
                       L.ForRange(i, 0, gdim, body=
                                  L.AssignSub(Xk[j], Kf[j,i]*(xk[i] - xgoal[i])))),
            ]

        newton_body = part1 + part2a  # Use if |dX| is needed
        #newton_body = part1 + part2b # Use for fixed iteration number

        # Loop until convergence
        num_iter = degree # TODO: Check if this is a good convergence criteria, add break if |dX| < epsilon?
        newton_loop = L.ForRange(k, 0, num_iter, body=newton_body)

        # X[ip] = Xk
        newton_finish = L.ForRange(j, 0, tdim, body=L.Assign(X[ip][j], Xk[j]))

        # Carry out newton loop for each point
        code = L.ForRange(ip, 0, num_points,
                          body=L.StatementList([newton_init, newton_loop, newton_finish]))
        return code

    def compute_jacobians(self, L, ir): # FIXME: Finish implementation of this, then mirror solution in compute_physical_coordinates
        # FIXME: Get data for scalar coordinate subelement:
        num_scalar_dofs = 3 # ir["num_scalar_coordinate_element_dofs"]
        scalar_coordinate_element_classname = "fixmecec" # ir["scalar_coordinate_element_classname"]

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
        # FIXME: Double check block structure of coordinate dofs
        coordinate_dofs = L.FlattenedArray(L.Symbol("coordinate_dofs"), dims=(num_scalar_dofs, gdim))

        # Output geometry
        J = L.FlattenedArray(L.Symbol("J"), dims=(num_points, gdim, tdim))

        # Input geometry
        X = L.FlattenedArray(L.Symbol("X"), dims=(num_points, tdim))

        # Declare basis derivatives table
        dphi = L.Symbol("dphi")
        dphi_dims = (tdim, num_scalar_dofs) # FIXME: Array layout to match eval_ref_bas_deriv
        dphi_decl = L.ArrayDecl("double", dphi, sizes=dphi_dims)

        # Computing table one point at a time instead of using
        # num_points will allow skipping dynamic allocation
        one_point = 1

        # Define scalar finite element instance (stateless, so placing this on the stack is free)
        # FIXME: To do this we'll need to #include the element header. Find a solution with dijitso!!
        #        When that's fixed, we have a solution for custom integrals as well.
        define_element = "%s element;" % (scalar_coordinate_element_classname,)
        func = "element.evaluate_reference_basis_derivatives" # FIXME: Use correct function to compute basis derivatives here

        # Compute basis derivatives table
        compute_dphi = L.Call(func, (dphi, one_point, L.AddressOf(X[ip, 0]))) # FIXME: eval_ref_bas_deriv signature

        # Make table more accessible with dimensions
        dphi = L.FlattenedArray(dphi, dims=dphi_dims)

        # Assign to J[ip][i][j] for each component i,j
        J_loop = L.AssignAdd(J[ip, i, j], coordinate_dofs[d, i]*dphi[j, d]) # FIXME: Array layout of dphi to match eval_ref_bas_deriv
        J_loop = L.ForRange(d, 0, num_scalar_dofs, body=J_loop)
        J_loop = L.ForRange(j, 0, tdim, body=J_loop)
        J_loop = L.ForRange(i, 0, gdim, body=J_loop)

        # Carry out computation of dphi and J accumulation for each point
        point_body = L.StatementList([compute_dphi, J_loop])
        point_loop = L.ForRange(ip, 0, num_points, body=point_body)

        body = L.StatementList([define_element, dphi_decl, point_loop])
        return body

    def compute_jacobian_determinants(self, L, ir): # Looks good on inspection. TODO: Test!
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

        # Assign det expression to detJ # TODO: Call Eigen instead?
        if gdim == tdim:
            body = L.Assign(detJ, det_nn(J[ip], gdim))
        elif tdim == 1:
            body = L.Assign(detJ, cell_orientation*pdet_m1(A, gdim))
        #elif tdim == 2 and gdim == 3:
        #    body = L.Assign(detJ, cell_orientation*pdet_32(A)) # Possible optimization not implemented here
        else:
            JTJ = L.Symbol("JTJ")
            body = L.StatementList([
                generate_compute_ATA(L, JTJ, J[ip], gdim, tdim),
                L.Assign(detJ, cell_orientation*L.Call("sqrt", det_nn(JTJ, tdim))),
                ])

        # Carry out for all points
        loop = L.ForRange(ip, 0, num_points, body=body)

        code = loop #[defines, loop]
        return code

    def compute_jacobian_inverses(self, L, ir): # FIXME: Not working yet. Call ufc_geometry function or eigen?
        # Dimensions
        gdim = ir["geometric_dimension"]
        tdim = ir["topological_dimension"]
        num_points = L.Symbol("num_points")

        # Loop indices
        ip = L.Symbol("ip")
        i = L.Symbol("i")
        j = L.Symbol("j")

        # Input cell data
        coordinate_dofs = L.Symbol("coordinate_dofs")
        cell_orientation = L.Symbol("cell_orientation")

        # Output geometry
        K = L.FlattenedArray(L.Symbol("K"), dims=(num_points, tdim, gdim))

        # Input geometry
        J = L.FlattenedArray(L.Symbol("J"), dims=(num_points, gdim, tdim))
        detJ = L.Symbol("detJ")

        # Assign to K[j][i] for each component j,i
        # FIXME: compute K[ip] from  J[ip], detJ[ip]
        body = L.Assign(K[ip][j][i], 0.0) # FIXME: Call Eigen?
        body = L.ForRange(i, 0, gdim, body=body)
        body = L.ForRange(j, 0, tdim, body=body)

        # Carry out for all points
        return L.ForRange(ip, 0, num_points, body=body)

    def compute_geometry(self, L, ir): # Output looks good on inspection. TODO: Test!
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
            L.Call("compute_physical_coordinates", (x, num_points, X, coordinate_dofs, cell_orientation)),
            L.Call("compute_jacobians", (J, num_points, X, coordinate_dofs, cell_orientation)),
            L.Call("compute_jacobian_determinants", (detJ, num_points, J)),
            L.Call("compute_jacobian_inverses", (K, num_points, J, detJ)),
            ]
        return L.StatementList(code)
