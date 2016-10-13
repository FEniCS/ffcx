# -*- coding: utf-8 -*-
"""Work in progress translation of FFC evaluatebasis code to uflacs CNodes format."""

from six import string_types
from ffc.log import error

"""
TODO: Add these to ufc::finite_element:

    /// NEW: Evaluate all reference basis function values at given points X in reference cell
    /// Unflattened shape of values is [num_points][num_dofs][reference_value_size]
    /// Unflattened shape of X is [num_points][tdim]
    virtual void evaluate_reference_basis_values(double * reference_values,
                                                 std::size_t num_points,
                                                 const double * X) const = 0;

    /// NEW: Evaluate all reference basis function derivatives up to order at given points X in reference cell
    /// Unflattened shape of values is [num_points][num_derivatives][num_dofs][reference_value_size]
    /// Unflattened shape of X is [num_points][tdim]
    virtual void evaluate_reference_basis_derivatives(double * reference_values,
                                                      std::size_t num_points,
                                                      std::size_t order,
                                                      const double * X) const = 0;

"""

''' /// FUTURE SPLIT IMPLEMENTATION OF EVALUATE_BASIS:
    /// Evaluate basis function i at given point x in cell
    virtual void evaluate_basis(std::size_t i,
                                double* values,
                                const double* x,
                                const double* coordinate_dofs,
                                int cell_orientation) const;
    /// Evaluate all basis functions at given point x in cell
    virtual void evaluate_basis_all(double* values,
                                    const double* x,
                                    const double* coordinate_dofs,
                                    int cell_orientation) const
    ... and derivatives
    {
      const std::size_t gdim = 3;
      const std::size_t tdim = 2;
      const std::size_t num_points = 1;

      // domain::
      double X[num_points*tdim]; // X[i] -> X[ip*tdim + i]
      compute_reference_coordinates(X, num_points, x, coordinate_dofs, cell_orientation);

      // domain::
      double J[num_points*gdim*tdim]; // J[i,j] -> J[ip*gdim*tdim + i*tdim + j]
      compute_jacobians(J, num_points, X, coordinate_dofs, cell_orientation);

      // domain::
      double detJ[num_points]; // detJ -> detJ[ip]
      compute_jacobian_determinants(detJ, num_points, J);

      // domain::
      double K[num_points*tdim*gdim]; // K[i,j] -> K[ip*tdim*gdim + i*gdim + j]
      compute_jacobian_inverses(K, num_points, J, detJ);

      // domain:: (inverse of compute_reference_coordinates)
      //double x[num_points*gdim]; // x[i] -> x[ip*gdim + i]
      //compute_physical_coordinates(x, num_points, X, K, coordinate_dofs, cell_orientation);

      // domain:: (combining the above)
      //compute_geometry(x, J, detJ, K, num_points, X, coordinate_dofs, cell_orientation);

      // phi[ip*ndofs*rvs + idof*rvs + jcomp]
      double reference_basis_values[num_points*num_dofs*reference_value_size];
      compute_reference_basis(reference_basis_values, num_points, X);

      // phi[ip*nder*ndofs*rvs + iderivative*ndofs*rvs + idof*rvs + jcomp]
      double reference_basis_derivatives[num_points*num_derivatives*num_dofs*reference_value_size];
      compute_reference_basis_derivatives(reference_basis_derivatives, derivative_order, num_points, X);

      double physical_basis_values[num_points*num_dofs*value_size]; // phi -> phi[ip*ndofs*pvs + idof*pvs + icomp]
      compute_physical_basis[_derivatives](physical_basis_values, num_points, reference_basis_values, J, detJ, K);
    }
'''

import math

def generate_evaluate_reference_basis(L, data):
    """Generate code to evaluate element basisfunctions at an arbitrary point on the reference element.

    The value(s) of the basisfunction is/are computed as in FIAT as
    the dot product of the coefficients (computed at compile time)
    and basisvalues which are dependent on the coordinate and thus
    have to be computed at run time.

    The function should work for all elements supported by FIAT, but
    it remains untested for tensor valued elements.

    This code is adapted from code in FFC which computed the basis
    from physical coordinates, and also to use UFLACS utilities.

    The FFC code has a comment "From FIAT_NEW.polynomial_set.tabulate()".
    """
    # Cutoff for feature to disable generation of this code (consider removing after benchmarking final result)
    if isinstance(data, string_types):
        return L.Throw("evaluate_reference_basis: %s" % data)

    # Get some known dimensions
    element_cellname = data["cellname"]
    #gdim = data["geometric_dimension"]
    tdim = data["topological_dimension"]
    reference_value_size = data["reference_value_size"]
    num_dofs = len(data["dofs_data"])

    # Input geometry
    num_points = L.Symbol("num_points")
    X = L.Symbol("X")

    # Output values
    reference_values = L.Symbol("reference_values")
    ref_values = L.FlattenedArray(reference_values, dims=(num_points, num_dofs, reference_value_size))

    # NB! This symbol refers to the FIAT reference coordinate!
    Y = L.Symbol("Y")

    # Loop indices
    ip = L.Symbol("ip")
    k = L.Symbol("k")
    c = L.Symbol("c")
    r = L.Symbol("r")


    # Generate code with static tables of expansion coefficients
    tables_code = []
    coefficients_for_dof = []
    for idof, dof_data in enumerate(data["dofs_data"]):
        num_components = dof_data["num_components"]
        num_members = dof_data["num_expansion_members"]
        fiat_coefficients = dof_data["coeffs"]

        # TODO: Check if any fiat_coefficients tables in expansion_coefficients are equal and reuse instead of declaring new.

        # Create separate variable name for coefficients table for each dof
        coefficients = L.Symbol("coefficients%d" % idof)

        # Create static table with expansion coefficients computed by FIAT compile time.
        tables_code += [L.ArrayDecl("static const double", coefficients,
                                    (num_components, num_members), values=fiat_coefficients)]

        # Store symbol reference for this dof
        coefficients_for_dof.append(coefficients)


    # Reset values[:] to 0
    reset_values_code = [
        L.ForRange(k, 0, num_points*num_dofs*reference_value_size, body=
            L.Assign(reference_values[k], 0.0))
        ]


    # Mapping from UFC reference cell coordinate X to FIAT reference cell coordinate Y
    fiat_coordinate_mapping = L.ArrayDecl("const double", Y, (tdim,),
                                          values=[2.0*X[ip*tdim + jj]-1.0 for jj in range(tdim)])


    # Generate code to compute tables of basisvalues
    basisvalues_code = []
    basisvalues_for_degree = {}
    for idof, dof_data in enumerate(data["dofs_data"]):
        embedded_degree = dof_data["embedded_degree"]
        num_members = dof_data["num_expansion_members"]

        if embedded_degree not in basisvalues_for_degree:
            basisvalues = L.Symbol("basisvalues%d" % embedded_degree)
            bfcode = _generate_compute_basisvalues(L, basisvalues, Y, element_cellname, embedded_degree, num_members)
            basisvalues_code += [L.StatementList(bfcode)]

            # Store symbol reference for this degree
            basisvalues_for_degree[embedded_degree] = basisvalues


    # Accumulate products of basisvalues and coefficients into values
    accumulation_code = []
    for idof, dof_data in enumerate(data["dofs_data"]):
        embedded_degree = dof_data["embedded_degree"]
        num_components = dof_data["num_components"]
        num_members = dof_data["num_expansion_members"]

        # In ffc representation, this is extracted per dof (but will coincide for some dofs of piola mapped elements):
        reference_offset = dof_data["reference_offset"]

        # Select the right basisvalues for this dof
        basisvalues = basisvalues_for_degree[embedded_degree]

        # Select the right coefficients for this dof
        coefficients = coefficients_for_dof[idof]

        # Generate basis accumulation loop
        if num_components > 1:
            accumulation_code += [
                L.ForRange(c, 0, num_components, body=
                    L.ForRange(r, 0, num_members, body=
                        L.AssignAdd(ref_values[ip, idof, reference_offset + c],
                                    coefficients[c,r] * basisvalues[r])))
                ]
        elif num_members > 1:
            accumulation_code += [
                L.ForRange(r, 0, num_members, body=
                    L.AssignAdd(ref_values[ip, idof, reference_offset], coefficients[0, r] * basisvalues[r]))
                ]
        else:
            accumulation_code += [
                L.AssignAdd(ref_values[ip, idof, reference_offset], coefficients[0, 0] * basisvalues[0])
                ]

        # FIXME: Move this mapping to its own ufc function e.g. finite_element::apply_element_mapping(reference_values, J, K)
        #code += _generate_apply_mapping_to_computed_values(L, dof_data) # Only works for affine (no-op)

    # Stitch it all together
    code = L.StatementList(
        tables_code +
        reset_values_code +
        [L.ForRange(ip, 0, num_points, body=L.StatementList([
            L.Comment("Map from UFC reference coordinate X to FIAT reference coordinate Y"),
            fiat_coordinate_mapping,
            L.Comment("Compute basisvalues for each relevant embedded degree"),
            basisvalues_code,
            L.Comment("Accumulate products of coefficients and basisvalues"),
            accumulation_code,
            ]))])
    return code


def _generate_compute_basisvalues(L, basisvalues, Y, element_cellname, embedded_degree, num_members):
    """From FIAT_NEW.expansions."""

    # Branch off to cell specific implementations
    if element_cellname == "interval":
        code = _generate_compute_interval_basisvalues(L, basisvalues, Y, embedded_degree, num_members)
    elif element_cellname == "triangle":
        code = _generate_compute_triangle_basisvalues(L, basisvalues, Y, embedded_degree, num_members)
    elif element_cellname == "tetrahedron":
        code = _generate_compute_tetrahedron_basisvalues(L, basisvalues, Y, embedded_degree, num_members)
    else:
        error()

    return code

def _jrc(a, b, n):
    an = float((2*n+1+a+b) * (2*n+2+a+b))   / float(2*(n+1) * (n+1+a+b))
    bn = float((a*a-b*b) * (2*n+1+a+b))     / float(2*(n+1) * (2*n+a+b) * (n+1+a+b))
    cn = float((n+a) * (n+b) * (2*n+2+a+b)) / float((n+1) * (n+1+a+b) * (2*n+a+b))
    return (an, bn, cn)

def _generate_compute_interval_basisvalues(L, basisvalues, Y, embedded_degree, num_members):
    # FIAT_NEW.expansions.LineExpansionSet.

    # Create zero-initialized array for with basisvalues
    code = [L.ArrayDecl("double", basisvalues, (num_members,), values=0)]

    # FIAT_NEW.jacobi.eval_jacobi_batch(a,b,n,xs)
    # for ii in range(result.shape[1]):
    #    result[0,ii] = 1.0 + xs[ii,0] - xs[ii,0]
    # The initial value basisvalues[0] is always 1.0
    code += [L.Assign(basisvalues[0], 1.0)]

    if embedded_degree > 0:
        # FIAT_NEW.jacobi.eval_jacobi_batch(a,b,n,xs).
        # result[1,:] = 0.5 * ( a - b + ( a + b + 2.0 ) * xsnew )
        # The initial value basisvalues[1] is always Y[0]
        code += [L.Assign(basisvalues[1], Y[0])]

    # Only active if embedded_degree > 1.
    for r in range(2, embedded_degree+1):
        # FIAT_NEW.jacobi.eval_jacobi_batch(a,b,n,xs).
        # apb = a + b (equal to 0 because of function arguments)
        # for k in range(2,n+1):
        #    a1 = 2.0 * k * ( k + apb ) * ( 2.0 * k + apb - 2.0 )
        #    a2 = ( 2.0 * k + apb - 1.0 ) * ( a * a - b * b )
        #    a3 = ( 2.0 * k + apb - 2.0 )  \
        #        * ( 2.0 * k + apb - 1.0 ) \
        #        * ( 2.0 * k + apb )
        #    a4 = 2.0 * ( k + a - 1.0 ) * ( k + b - 1.0 ) \
        #        * ( 2.0 * k + apb )
        #    a2 = a2 / a1
        #    a3 = a3 / a1
        #    a4 = a4 / a1
        #    result[k,:] = ( a2 + a3 * xsnew ) * result[k-1,:] \
        #        - a4 * result[k-2,:]

        # The below implements the above (with a = b = apb = 0)
        a1 = float(2*r*r*(2*r - 2))
        a3 = ((2*r - 2)*(2*r - 1)*(2*r)) / a1
        a4 = (2*(r - 1)*(r - 1)*(2*r)) / a1
        value = (Y[0] * (a3 * basisvalues[r-1])) - a4*basisvalues[r-2]
        code += [L.Assign(basisvalues[r], value)]

    # FIAT_NEW code
    # results = numpy.zeros( ( n+1 , len(pts) ) , type( pts[0][0] ) )
    # for k in range( n + 1 ):
    #    results[k,:] = psitilde_as[k,:] * math.sqrt( k + 0.5 )

    # Scale values
    p = L.Symbol("p")
    code += [L.ForRange(p, 0, embedded_degree + 1,
                        body=L.AssignMul(basisvalues[p], L.Sqrt(0.5 + p)))]
    return code

def _generate_compute_triangle_basisvalues(L, basisvalues, Y, embedded_degree, num_members):
    # FIAT_NEW.expansions.TriangleExpansionSet.

    def _idx2d(p, q):
        return (p + q)*(p + q + 1)//2 + q

    # Create zero-initialized array for with basisvalues
    code = [L.ArrayDecl("double", basisvalues, (num_members,), values=0)]

    # Compute helper factors
    # FIAT_NEW code
    # f1 = (1.0+2*x+y)/2.0
    # f2 = (1.0 - y) / 2.0
    # f3 = f2**2

    # The initial value basisvalues 0 is always 1.0.
    # FIAT_NEW code
    # for ii in range( results.shape[1] ):
    #    results[0,ii] = 1.0 + apts[ii,0]-apts[ii,0]+apts[ii,1]-apts[ii,1]
    code += [L.Assign(basisvalues[0], 1.0)]

    # Only continue if the embedded degree is larger than zero.
    if embedded_degree == 0:
        return code

    # The initial value of basisfunction 1 is equal to f1.
    # FIAT_NEW code
    # results[idx(1,0),:] = f1
    f1 = L.Symbol("tmp1_%d" % embedded_degree)
    code += [L.VariableDecl("const double", f1, (1.0 + 2.0*Y[0] + Y[1]) / 2.0)]
    code += [L.Assign(basisvalues[1], f1)]

    # NOTE: KBO: The order of the loops is VERY IMPORTANT!!

    # FIAT_NEW code (loop 1 in FIAT)
    # for p in range(1,n):
    #    a = (2.0*p+1)/(1.0+p)
    #    b = p / (p+1.0)
    #    results[idx(p+1,0)] = a * f1 * results[idx(p,0),:] \
    #        - b * f3 *results[idx(p-1,0),:]
    # Only active if embedded_degree > 1.
    if embedded_degree > 1:
        f2 = L.Symbol("tmp2_%d" % embedded_degree)
        f3 = L.Symbol("tmp3_%d" % embedded_degree)
        code += [L.VariableDecl("const double", f2, (1.0 - Y[1]) / 2.0)]
        code += [L.VariableDecl("const double", f3, f2*f2)]
    for r in range(1, embedded_degree):
        rr = _idx2d((r + 1), 0)
        ss = _idx2d(r, 0)
        tt = _idx2d((r - 1), 0)
        A = (2*r + 1.0)/(r + 1)
        B = r/(1.0 + r)
        v1 = A * f1 * basisvalues[ss]
        v2 = B * f3 * basisvalues[tt]
        value = v1 - v2
        code += [L.Assign(basisvalues[rr], value)]

    # FIAT_NEW code (loop 2 in FIAT).
    # for p in range(n):
    #    results[idx(p,1),:] = 0.5 * (1+2.0*p+(3.0+2.0*p)*y) \
    #        * results[idx(p,0)]
    for r in range(0, embedded_degree):
        # (p+q)*(p+q+1)//2 + q
        rr = _idx2d(r, 1)
        ss = _idx2d(r, 0)
        A = 0.5*(1 + 2*r)
        B = 0.5*(3 + 2*r)
        C = A + (B * Y[1])
        value = C * basisvalues[ss]
        code += [L.Assign(basisvalues[rr], value)]

    # Only active if embedded_degree > 1.
    # FIAT_NEW code (loop 3 in FIAT).
    # for p in range(n-1):
    #    for q in range(1,n-p):
    #        (a1,a2,a3) = jrc(2*p+1,0,q)
    #        results[idx(p,q+1),:] \
    #            = ( a1 * y + a2 ) * results[idx(p,q)] \
    #            - a3 * results[idx(p,q-1)]
    # Only active if embedded_degree > 1.
    for r in range(0, embedded_degree - 1):
        for s in range(1, embedded_degree - r):
            rr = _idx2d(r, (s + 1))
            ss = _idx2d(r, s)
            tt = _idx2d(r, s - 1)
            A, B, C = _jrc(2*r + 1, 0, s)
            value = (B + A*Y[1])*basisvalues[ss] - C*basisvalues[tt]
            code += [L.Assign(basisvalues[rr], value)]

    # FIAT_NEW code (loop 4 in FIAT).
    # for p in range(n+1):
    #    for q in range(n-p+1):
    #        results[idx(p,q),:] *= math.sqrt((p+0.5)*(p+q+1.0))
    for r in range(0, embedded_degree + 1):
        for s in range(0, embedded_degree + 1 - r):
            rr = _idx2d(r, s)
            A = (r + 0.5)*(r + s + 1)
            code += [L.AssignMul(basisvalues[rr], math.sqrt(A))]

    return code

def _generate_compute_tetrahedron_basisvalues(L, basisvalues, Y, embedded_degree, num_members):
    # FIAT_NEW.expansions.TetrahedronExpansionSet.

    # FIAT_NEW code (compute index function) TetrahedronExpansionSet.
    # def idx(p,q,r):
    #     return (p+q+r)*(p+q+r+1)*(p+q+r+2)//6 + (q+r)*(q+r+1)//2 + r
    def _idx3d(p, q, r):
        return (p+q+r)*(p+q+r+1)*(p+q+r+2)//6 + (q+r)*(q+r+1)//2 + r

    # Compute helper factors.
    # FIAT_NEW code
    # factor1 = 0.5 * ( 2.0 + 2.0*x + y + z )
    # factor2 = (0.5*(y+z))**2
    # factor3 = 0.5 * ( 1 + 2.0 * y + z )
    # factor4 = 0.5 * ( 1 - z )
    # factor5 = factor4 ** 2

    # Create zero-initialized array for with basisvalues
    code = [L.ArrayDecl("double", basisvalues, (num_members,), values=0)]

    # The initial value basisvalues 0 is always 1.0.
    # FIAT_NEW code
    # for ii in range( results.shape[1] ):
    #    results[0,ii] = 1.0 + apts[ii,0]-apts[ii,0]+apts[ii,1]-apts[ii,1]
    code += [L.Assign(basisvalues[0], 1.0)]

    # Only continue if the embedded degree is larger than zero.
    if embedded_degree == 0:
        return code

    # The initial value of basisfunction 1 is equal to f1.
    # FIAT_NEW code
    # results[idx(1,0),:] = f1
    f1 = L.Symbol("tmp1_%d" % embedded_degree)
    code += [L.VariableDecl("const double", f1, 0.5*(2.0 + 2.0*Y[0] + Y[1] + Y[2]))]
    code += [L.Assign(basisvalues[1], f1)]

    # NOTE: KBO: The order of the loops is VERY IMPORTANT!!

    # FIAT_NEW code (loop 1 in FIAT).
    # for p in range(1,n):
    #    a1 = ( 2.0 * p + 1.0 ) / ( p + 1.0 )
    #    a2 = p / (p + 1.0)
    #    results[idx(p+1,0,0)] = a1 * factor1 * results[idx(p,0,0)] \
    #        -a2 * factor2 * results[ idx(p-1,0,0) ]
    # Only active if embedded_degree > 1.
    if embedded_degree > 1:
        f2 = L.Symbol("tmp2_%d" % embedded_degree)
        code += [L.VariableDecl("const double", f2, 0.25*(Y[1] + Y[2])*(Y[1] + Y[2]))]
    for r in range(1, embedded_degree):
        rr = _idx3d((r + 1), 0, 0)
        ss = _idx3d(r, 0, 0)
        tt = _idx3d((r - 1), 0, 0)
        A = (2*r + 1.0)/(r + 1)
        B = r/(r + 1.0)
        value = (A * f1 * basisvalues[ss]) - (B * f2 * basisvalues[tt])
        code += [L.Assign(basisvalues[rr], value)]

    # FIAT_NEW code (loop 2 in FIAT).
    # q = 1
    # for p in range(0,n):
    #    results[idx(p,1,0)] = results[idx(p,0,0)] \
    #        * ( p * (1.0 + y) + ( 2.0 + 3.0 * y + z ) / 2 )
    for r in range(0, embedded_degree):
        rr = _idx3d(r, 1, 0)
        ss = _idx3d(r, 0, 0)
        term0 = 0.5 * (2.0 + 3.0*Y[1] + Y[2])
        term1 = float(r) * (1.0 + Y[1])
        value = (term0 + term1) * basisvalues[ss]
        code += [L.Assign(basisvalues[rr], value)]

    # FIAT_NEW code (loop 3 in FIAT).
    # for p in range(0,n-1):
    #    for q in range(1,n-p):
    #        (aq,bq,cq) = jrc(2*p+1,0,q)
    #        qmcoeff = aq * factor3 + bq * factor4
    #        qm1coeff = cq * factor5
    #        results[idx(p,q+1,0)] = qmcoeff * results[idx(p,q,0)] \
    #            - qm1coeff * results[idx(p,q-1,0)]
    # Only active if embedded_degree > 1.
    if embedded_degree > 1:
        f3 = L.Symbol("tmp3_%d" % embedded_degree)
        f4 = L.Symbol("tmp4_%d" % embedded_degree)
        f5 = L.Symbol("tmp5_%d" % embedded_degree)
        code += [L.VariableDecl("const double", f3, 0.5*(1.0 + 2.0*Y[1] + Y[2]))]
        code += [L.VariableDecl("const double", f4, 0.5*(1.0 - Y[2]))]
        code += [L.VariableDecl("const double", f5, f4*f4)]
    for r in range(0, embedded_degree - 1):
        for s in range(1, embedded_degree - r):
            rr = _idx3d(r, (s + 1), 0)
            ss = _idx3d(r, s, 0)
            tt = _idx3d(r, s - 1, 0)
            (A, B, C) = _jrc(2*r + 1, 0, s)
            term0 = ((A*f3) + (B*f4)) * basisvalues[ss]
            term1 = C * f5 * basisvalues[tt]
            value = term0 - term1
            code += [L.Assign(basisvalues[rr], value)]

    # FIAT_NEW code (loop 4 in FIAT).
    # now handle r=1
    # for p in range(n):
    #    for q in range(n-p):
    #        results[idx(p,q,1)] = results[idx(p,q,0)] \
    #            * ( 1.0 + p + q + ( 2.0 + q + p ) * z )
    for r in range(0, embedded_degree):
        for s in range(0, embedded_degree - r):
            rr = _idx3d(r, s, 1)
            ss = _idx3d(r, s, 0)
            A = (float(2 + r + s) * Y[2]) + float(1 + r + s)
            value = A * basisvalues[ss]
            code += [L.Assign(basisvalues[rr], value)]

    # FIAT_NEW code (loop 5 in FIAT).
    # general r by recurrence
    # for p in range(n-1):
    #     for q in range(0,n-p-1):
    #         for r in range(1,n-p-q):
    #             ar,br,cr = jrc(2*p+2*q+2,0,r)
    #             results[idx(p,q,r+1)] = \
    #                         (ar * z + br) * results[idx(p,q,r) ] \
    #                         - cr * results[idx(p,q,r-1) ]
    # Only active if embedded_degree > 1.
    for r in range(0, embedded_degree - 1):
        for s in range(0, embedded_degree - r - 1):
            for t in range(1, embedded_degree - r - s):
                rr = _idx3d(r, s, ( t + 1))
                ss = _idx3d(r, s, t)
                tt = _idx3d(r, s, t - 1)
                (A, B, C) = _jrc(2*r + 2*s + 2, 0, t)
                az_b = B + A*Y[2]
                value = (az_b*basisvalues[ss]) - (C*basisvalues[tt])
                code += [L.Assign(basisvalues[rr], value)]

    # FIAT_NEW code (loop 6 in FIAT).
    # for p in range(n+1):
    #    for q in range(n-p+1):
    #        for r in range(n-p-q+1):
    #            results[idx(p,q,r)] *= math.sqrt((p+0.5)*(p+q+1.0)*(p+q+r+1.5))
    for r in range(0, embedded_degree + 1):
        for s in range(0, embedded_degree - r + 1):
            for t in range(0, embedded_degree - r - s + 1):
                rr = _idx3d(r, s, t)
                A = (r + 0.5)*(r + s + 1)*(r + s + t + 1.5)
                code += [L.AssignMul(basisvalues[rr], math.sqrt(A))]

    return code


def _generate_apply_mapping_to_computed_values(L, dof_data):
    mapping = dof_data["mapping"]
    num_components = dof_data["num_components"]
    reference_offset = dof_data["reference_offset"]
    physical_offset = dof_data["physical_offset"]

    physical_values[num_points][num_dofs][physical_value_size]
    reference_values[num_points][num_dofs][reference_value_size]

    code = []

    # FIXME: Define values numbering
    if mapping == "affine":
        # Just copy values
        if num_components == 1:
            code += [
                L.ForRange(ip, 0, num_points, body=
                    L.Assign(physical_values[ip][physical_offset],
                             reference_values[ip][reference_offset])),
                ]
        else:
            code += [
                L.ForRange(ip, 0, num_points, body=
                    L.ForRange(k, 0, num_components, body=
                        L.Assign(physical_values[ip][physical_offset + k],
                                 reference_values[ip][reference_offset + k]))),
                ]
        return code

    else:
        fixme

    return code

def __ffc_implementation_of__generate_apply_mapping_to_computed_values(L):
    # Apply transformation if applicable.
    mapping = dof_data["mapping"]
    num_components = dof_data["num_components"]
    reference_offset = dof_data["reference_offset"]
    physical_offset = dof_data["physical_offset"]

    if mapping == "affine":
        pass

    elif mapping == "contravariant piola":
        code += ["", f_comment("Using contravariant Piola transform to map values back to the physical element")]

        # Get temporary values before mapping.
        code += [f_const_float(f_tmp_ref(i), f_component(f_values, i + offset))
                 for i in range(num_components)]

        # Create names for inner product.
        basis_col = [f_tmp_ref(j) for j in range(tdim)]
        for i in range(gdim):
            # Create Jacobian.
            jacobian_row = [f_trans("J", i, j, gdim, tdim, None) for j in range(tdim)]
            # Create inner product and multiply by inverse of Jacobian.
            inner = f_group(f_inner(jacobian_row, basis_col)) # sum_j J[i,j], values[offset + j])
            value = f_mul([f_inv(f_detJ(None)), inner])
            name = f_component(f_values, i + offset)
            # FIXME: This is writing values[offset+:] = M[:,:] * values[offset+:],
            #        i.e. offset must be physical (unless there's a bug).
            #        We want to use reference offset for evaluate_reference_basis,
            #        and to make this mapping read from one array using reference offset
            #        and write to another array using physical offset.
            code += [f_assign(name, value)]

    elif mapping == "covariant piola":
        code += ["", f_comment("Using covariant Piola transform to map values back to the physical element")]
        # Get temporary values before mapping.
        code += [f_const_float(f_tmp_ref(i), f_component(f_values, i + offset))
                 for i in range(num_components)]
        # Create names for inner product.
        tdim = data["topological_dimension"]
        gdim = data["geometric_dimension"]
        basis_col = [f_tmp_ref(j) for j in range(tdim)]
        for i in range(gdim):
            # Create inverse of Jacobian.
            inv_jacobian_column = [f_trans("JINV", j, i, tdim, gdim, None) for j in range(tdim)]

            # Create inner product of basis values and inverse of Jacobian.
            value = f_group(f_inner(inv_jacobian_column, basis_col))
            name = f_component(f_values, i + offset)
            code += [f_assign(name, value)]

    elif mapping == "double covariant piola":
        code += ["", f_comment("Using double covariant Piola transform to map values back to the physical element")]
        # Get temporary values before mapping.
        code += [f_const_float(f_tmp_ref(i), f_component(f_values, i + offset))
                 for i in range(num_components)]
        # Create names for inner product.
        tdim = data["topological_dimension"]
        gdim = data["geometric_dimension"]
        basis_col = [f_tmp_ref(j) for j in range(num_components)]
        for p in range(num_components):
            # unflatten the indices
            i = p // tdim
            l = p % tdim
            # g_il = K_ji G_jk K_kl
            value = f_group(f_inner(
                [f_inner([f_trans("JINV", j, i, tdim, gdim, None)
                          for j in range(tdim)],
                         [basis_col[j * tdim + k] for j in range(tdim)])
                 for k in range(tdim)],
                [f_trans("JINV", k, l, tdim, gdim, None)
                 for k in range(tdim)]))
            name = f_component(f_values, p + offset)
            code += [f_assign(name, value)]

    elif mapping == "double contravariant piola":
        code += ["", f_comment("Pullback of a matrix-valued funciton as contravariant 2-tensor mapping values back to the physical element")]
        # Get temporary values before mapping.
        code += [f_const_float(f_tmp_ref(i), f_component(f_values, i + offset))
                 for i in range(num_components)]
        # Create names for inner product.
        tdim = data["topological_dimension"]
        gdim = data["geometric_dimension"]
        basis_col = [f_tmp_ref(j) for j in range(num_components)]
        for p in range(num_components):
            # unflatten the indices
            i = p // tdim
            l = p % tdim
            # g_il = (detJ)^(-2) J_ij G_jk J_lk
            value = f_group(f_inner(
                [f_inner([f_trans("J", i, j, tdim, gdim, None)
                          for j in range(tdim)],
                         [basis_col[j * tdim + k] for j in range(tdim)])
                 for k in range(tdim)],
                [f_trans("J", l, k, tdim, gdim, None)
                 for k in range(tdim)]))
            value = f_mul([f_inv(f_detJ(None)), f_inv(f_detJ(None)), value])
            name = f_component(f_values, p + offset)
            code += [f_assign(name, value)]

    else:
        error("Unknown mapping: %s" % mapping)
