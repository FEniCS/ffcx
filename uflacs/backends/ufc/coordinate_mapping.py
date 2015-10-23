
from uflacs.backends.ufc.generator import ufc_generator

### Code generation utilities:

def flat_array(L, name, dims):
    return L.FlattenedArray(L.Symbol(name), dims=dims)

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

def __pdet_m1(L, A, m):
    # Special case 1xm for simpler expression
    i = L.Symbol("i")
    A2 = A[i,0]*A[i,0] # TODO: Translate to code
    return L.Call("sqrt", A2)

def __pdet_23(L, A):
    # Special case 2x3 for simpler expression
    i = L.Symbol("i")

    # TODO: Translate to code:
    c = cross_expr(A[:,0], A[:,1])
    c2 = c[i]*c[i]

    return L.Call("sqrt", c2)

def pdet_mn(A, m, n):
    """Compute the pseudo-determinant of A: sqrt(det(A.T*A))."""
    # TODO: This would be more reusable if it didn't make up variable names...
    # Build A^T*A matrix
    i = L.Symbol("i")
    j = L.Symbol("j")
    k = L.Symbol("k")
    ATA = L.ArrayDecl("double", "ATA", shape=(n, n), values=0)
    body = L.AssignAdd(ATA[i, j], A[k, i] * A[k, j])
    body = L.ForRange(k, 0, m, body=body)
    body = L.ForRange(j, 0, n, body=body)
    body = L.ForRange(i, 0, n, body=body)

    # Take determinant and square root
    return L.Call("sqrt", det_nn(ATA, n))

def pdet_expr(A, m, n):
    """Compute the pseudo-determinant of A: sqrt(det(A.T*A))."""
    if n == 1:
        return pdet_mn(A, m, n)
        #return pdet_m1(A, m)
    elif m == 3 and n == 2:
        return pdet_mn(A, m, n)
        #return pdet_32(A)
    else:
        return pdet_mn(A, m, n)

def det_expr(A, m, n):
    "Compute the (pseudo-)determinant of A."
    if m == n:
        return det_nn(A, m)
    else:
        return pdet_expr(A, m, n)


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

    def compute_physical_coordinates(self, L, ir):
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
        # Dimensions
        gdim = ir["geometric_dimension"]
        tdim = ir["topological_dimension"]
        num_points = L.Symbol("num_points")

        # Computing table one point at a time instead of using
        # num_points will allow skipping dynamic allocation
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
        X = flat_array(L, "X", (num_points, tdim))

        # Input geometry
        x = flat_array(L, "x", (num_points, gdim))

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

        # TODO: Use cell midpoint instead of 0.0.
        initial_X_value = 0.0

        # Intermediate arrays declared in point loop
        Xk = L.Symbol("Xk")
        Xk_decl = L.ArrayDecl("double", Xk, (tdim,), values=initial_X_value)
        xgoal = L.Symbol("xgoal")
        xgoal_decl = L.ArrayDecl("double", xgoal, (gdim,), values=[x[ip][iv] for iv in range(gdim)])

        # Intermediate arrays declared in newton loop
        dX = L.Symbol("dX")
        dX_decl = L.ArrayDecl("double", dX, (tdim,), values=0.0) # Initialized to zero here!

        # Intermediate arrays to hold results of compute_geometry call
        xk = L.Symbol("xk")
        J = L.Symbol("J")
        detJ = L.Symbol("detJ")
        K = L.Symbol("K")
        xk_decl = L.ArrayDecl("double", xk, (gdim,))
        J_decl = L.ArrayDecl("double", J, (gdim*tdim,))
        detJ_decl = L.ArrayDecl("double", detJ, (1,))
        K_decl = L.ArrayDecl("double", K, (tdim*gdim,))

        # Compute K = J^-1 for one point (J and detJ are only used as
        # intermediate storage inside compute_geometry, not used out here)
        comp_geom = L.Call("compute_geometry",
                           (xk, J, detJ, K, one_point, Xk, coordinate_dofs, cell_orientation))

        # Compute dX[j] = sum_i K_ji * (x(Xk)_i - x_i)
        K = L.FlattenedArray(K, dims=(tdim, gdim))
        comp_dX = L.ForRange(j, 0, tdim, body=
                             L.ForRange(i, 0, gdim, body=
                                        L.AssignAdd(dX[j], K[j,i]*(xk[i] - xgoal[i]))))

        # Update Xk -= dX
        update_Xk = L.ForRange(j, 0, tdim, body=L.AssignSub(Xk[j], dX[j]))

        newton_body = [
            xk_decl, J_decl, detJ_decl, K_decl,
            comp_geom,
            dX_decl,
            comp_dX,
            update_Xk
            ]

        # declare xgoal = x[ip]
        # declare Xk = initial value
        newton_init = [
            xgoal_decl,
            Xk_decl,
            ]

        # Loop until convergence
        # TODO: Test convergence criteria
        num_iter = 3
        newton_loop = L.ForRange(k, 0, num_iter, body=newton_body)

        # X[ip] = Xk
        newton_finish = [
            L.ForRange(j, 0, tdim, body=L.Assign(X[ip][j], Xk[j]))
            ]

        # Join pieces
        code = L.StatementList([newton_init, newton_loop, newton_finish])

        # Carry out for all points
        return L.ForRange(ip, 0, num_points, body=code)

    def compute_jacobians(self, L, ir):
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
        coordinate_dofs = flat_array(L, "coordinate_dofs", (num_scalar_dofs, gdim)) # FIXME: Correct block structure of dofs?
        cell_orientation = L.Symbol("cell_orientation") # need this?

        # Output geometry
        J = flat_array(L, "J", (num_points, gdim, tdim))

        # Input geometry
        X = flat_array(L, "X", (num_points, tdim))

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


    def compute_jacobian_determinants(self, L, ir):
        # Dimensions
        gdim = ir["geometric_dimension"]
        tdim = ir["topological_dimension"]
        num_points = L.Symbol("num_points")

        # Loop indices
        ip = L.Symbol("ip")
        i = L.Symbol("i")
        j = L.Symbol("j")

        # Output geometry
        detJ = L.Symbol("detJ")[ip]

        # Input geometry
        J = flat_array(L, "J", (num_points, gdim, tdim))[ip]

        # Assign to detJ
        body = L.Assign(detJ, det_expr(J, gdim, tdim)) # TODO: Call Eigen instead?

        # Carry out for all points
        loop = L.ForRange(ip, 0, num_points, body=body)

        code = loop #[defines, loop]
        return code

    def compute_jacobian_inverses(self, L, ir):
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
        K = flat_array(L, "K", (num_points, tdim, gdim))[ip]

        # Input geometry
        J = flat_array(L, "J", (num_points, gdim, tdim))[ip]
        detJ = L.Symbol("detJ")[ip]

        # Assign to K[j][i] for each component j,i
        body = L.Assign(K[j][i], 0.0) # FIXME: Call Eigen?
        body = L.ForRange(i, 0, gdim, body=body)
        body = L.ForRange(j, 0, tdim, body=body)

        # Carry out for all points
        return L.ForRange(ip, 0, num_points, body=body)

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
            L.Call("compute_physical_coordinates", (x, num_points, X, coordinate_dofs, cell_orientation)),
            L.Call("compute_jacobians", (J, num_points, X, coordinate_dofs, cell_orientation)),
            L.Call("compute_jacobian_determinants", (detJ, num_points, J)),
            L.Call("compute_jacobian_inverses", (K, num_points, J, detJ)),
            ]
        return L.StatementList(code)
