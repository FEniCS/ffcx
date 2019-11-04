import logging

import numpy

logger = logging.getLogger(__name__)


# Used for various indices and arrays in this file
index_type = "int64_t"


def generate_evaluate_reference_basis(L, data, parameters):
    """Evaluate basis functions on the reference element.

    Generate code to evaluate element basisfunctions at an arbitrary
    point on the reference element.

    The value(s) of the basisfunction is/are computed as in FIAT as
    the dot product of the coefficients (computed at compile time) and
    basisvalues which are dependent on the coordinate and thus have to
    be computed at run time.

    The function should work for all elements supported by FIAT, but
    it remains untested for tensor valued elements.

    This code is adapted from code in FFC which computed the basis
    from physical coordinates, and also to use UFLACS utilities.

    The FFC code has a comment "From FIAT_NEW.polynomial_set.tabulate()".

    """
    # Cutoff for feature to disable generation of this code (consider
    # removing after benchmarking final result)
    if isinstance(data, str):
        # Return error code
        return [L.Return(-1)]

    # Get some known dimensions
    element_cellname = data["cellname"]
    tdim = data["topological_dimension"]
    reference_value_size = data["reference_value_size"]
    num_dofs = len(data["dofs_data"])

    # Input geometry
    num_points = L.Symbol("num_points")
    X = L.Symbol("X")

    # Output values
    reference_values = L.Symbol("reference_values")
    ref_values = L.FlattenedArray(
        reference_values, dims=(num_points, num_dofs, reference_value_size))

    # Loop indices
    ip = L.Symbol("ip")
    k = L.Symbol("k")
    c = L.Symbol("c")
    r = L.Symbol("r")

    # Return statement (value indicates that function is implemented)
    ret = L.Return(0)

    # Generate code with static tables of expansion coefficients
    tables_code, coefficients_for_dof = generate_expansion_coefficients(
        L, data["dofs_data"])

    # Reset reference_values[:] to 0
    reset_values_code = [
        L.ForRange(
            k,
            0,
            num_points * num_dofs * reference_value_size,
            index_type=index_type,
            body=L.Assign(reference_values[k], 0.0))
    ]
    setup_code = tables_code + reset_values_code

    # Generate code to compute tables of basisvalues
    basisvalues_code, basisvalues_for_degree, need_fiat_coordinates = \
        generate_compute_basisvalues(
            L, data["dofs_data"], element_cellname, tdim, X, ip)

    # Accumulate products of basisvalues and coefficients into values
    accumulation_code = [
        L.Comment("Accumulate products of coefficients and basisvalues"),
    ]
    for idof, dof_data in enumerate(data["dofs_data"]):
        embedded_degree = dof_data["embedded_degree"]
        num_components = dof_data["num_components"]
        num_members = dof_data["num_expansion_members"]

        # In ffc representation, this is extracted per dof
        # (but will coincide for some dofs of piola mapped elements):
        reference_offset = dof_data["reference_offset"]

        # Select the right basisvalues for this dof
        basisvalues = basisvalues_for_degree[embedded_degree]

        # Select the right coefficients for this dof
        coefficients = coefficients_for_dof[idof]

        # Generate basis accumulation loop
        if num_components > 1:
            # Could just simplify by using this generic code
            # and dropping the below two special cases
            accumulation_code += [
                L.ForRange(
                    c,
                    0,
                    num_components,
                    index_type=index_type,
                    body=L.ForRange(
                        r,
                        0,
                        num_members,
                        index_type=index_type,
                        body=L.AssignAdd(
                            ref_values[ip, idof, reference_offset + c],
                            coefficients[c, r] * basisvalues[r])))
            ]
        elif num_members > 1:
            accumulation_code += [
                L.ForRange(
                    r,
                    0,
                    num_members,
                    index_type=index_type,
                    body=L.AssignAdd(ref_values[ip, idof, reference_offset],
                                     coefficients[0, r] * basisvalues[r]))
            ]
        else:
            accumulation_code += [
                L.AssignAdd(ref_values[ip, idof, reference_offset],
                            coefficients[0, 0] * basisvalues[0])
            ]

        # TODO: Move this mapping to its own ufc function
        # e.g. finite_element::apply_element_mapping(reference_values,
        # J, K)
        # code += _generate_apply_mapping_to_computed_values(L, dof_data) # Only works for affine (no-op)

    # Stitch it all together
    code = [
        setup_code,
        L.ForRange(
            ip,
            0,
            num_points,
            index_type=index_type,
            body=basisvalues_code + accumulation_code), ret
    ]
    return code


def generate_expansion_coefficients(L, dofs_data):
    # TODO: Use precision parameter to format coefficients to match
    # legacy implementation and make regression tests more robust

    all_tables = []
    tables_code = []
    coefficients_for_dof = []
    for idof, dof_data in enumerate(dofs_data):
        num_components = dof_data["num_components"]
        num_members = dof_data["num_expansion_members"]
        fiat_coefficients = dof_data["coeffs"]

        # Check if any fiat_coefficients tables in expansion_coefficients
        # are equal and reuse instead of declaring new.
        coefficients = None
        # NB: O(n^2) loop over tables
        for symbol, table in all_tables:
            if table.shape == fiat_coefficients.shape and numpy.allclose(
                    table, fiat_coefficients):
                coefficients = symbol
                break

        # Create separate variable name for coefficients table for each dof
        if coefficients is None:
            coefficients = L.Symbol("coefficients%d" % idof)
            all_tables.append((coefficients, fiat_coefficients))

            # Create static table with expansion coefficients computed by FIAT compile time.
            tables_code += [
                L.ArrayDecl(
                    "static const double",
                    coefficients, (num_components, num_members),
                    values=fiat_coefficients)
            ]

        # Store symbol reference for this dof
        coefficients_for_dof.append(coefficients)

    return tables_code, coefficients_for_dof


def generate_compute_basisvalues(L, dofs_data, element_cellname, tdim, X, ip):
    basisvalues_code = [
        L.Comment("Compute basisvalues for each relevant embedded degree"),
    ]
    basisvalues_for_degree = {}
    need_fiat_coordinates = False
    Y = L.Symbol("Y")
    for idof, dof_data in enumerate(dofs_data):
        embedded_degree = dof_data["embedded_degree"]

        if embedded_degree:
            need_fiat_coordinates = True

        if embedded_degree not in basisvalues_for_degree:
            num_members = dof_data["num_expansion_members"]

            basisvalues = L.Symbol("basisvalues%d" % embedded_degree)
            bfcode = _generate_compute_basisvalues(
                L, basisvalues, Y, element_cellname, embedded_degree,
                num_members)
            basisvalues_code += [L.StatementList(bfcode)]

            # Store symbol reference for this degree
            basisvalues_for_degree[embedded_degree] = basisvalues

    if need_fiat_coordinates:
        # Mapping from UFC reference cell coordinate X to FIAT reference cell coordinate Y
        fiat_coordinate_mapping = [
            L.Comment(
                "Map from UFC reference coordinate X to FIAT reference coordinate Y"
            ),
            L.ArrayDecl(
                "const double",
                Y, (tdim, ),
                values=[2.0 * X[ip * tdim + jj] - 1.0 for jj in range(tdim)]),
        ]
        basisvalues_code = fiat_coordinate_mapping + basisvalues_code

    return basisvalues_code, basisvalues_for_degree, need_fiat_coordinates


def _generate_compute_basisvalues(L, basisvalues, Y, element_cellname,
                                  embedded_degree, num_members):
    """From FIAT_NEW.expansions."""

    # Branch off to cell specific implementations
    if element_cellname == "interval":
        code = _generate_compute_interval_basisvalues(
            L, basisvalues, Y, embedded_degree, num_members)
    elif element_cellname == "triangle":
        code = _generate_compute_triangle_basisvalues(
            L, basisvalues, Y, embedded_degree, num_members)
    elif element_cellname == "tetrahedron":
        code = _generate_compute_tetrahedron_basisvalues(
            L, basisvalues, Y, embedded_degree, num_members)
    elif element_cellname == "quadrilateral":
        code = _generate_compute_quad_basisvalues(
            L, basisvalues, Y, embedded_degree, num_members)
    elif element_cellname == "hexahedron":
        code = _generate_compute_hex_basisvalues(
            L, basisvalues, Y, embedded_degree, num_members)
    else:
        raise RuntimeError("Not supported:" + element_cellname)

    return code


def _jrc(a, b, n):
    an = float((2 * n + 1 + a + b) * (2 * n + 2 + a + b)) / float(2 * (n + 1) * (n + 1 + a + b))
    bn = float((a * a - b * b) * (2 * n + 1 + a + b)) / float(2 * (n + 1) * (2 * n + a + b) * (n + 1 + a + b))
    cn = float((n + a) * (n + b) * (2 * n + 2 + a + b)) / float((n + 1) * (n + 1 + a + b) * (2 * n + a + b))
    return (an, bn, cn)


def _generate_compute_interval_basisvalues(L, basisvalues, Y, embedded_degree,
                                           num_members):
    # FIAT_NEW.expansions.LineExpansionSet.

    # Create zero-initialized array for with basisvalues
    code = [L.ArrayDecl("double", basisvalues, (num_members, ), values=0)]

    code += [L.Assign(basisvalues[0], 1.0)]

    if embedded_degree > 0:
        code += [L.Assign(basisvalues[1], Y[0])]

    # Only active if embedded_degree > 1.
    for r in range(2, embedded_degree + 1):
        a1 = float(2 * r * r * (2 * r - 2))
        a3 = ((2 * r - 2) * (2 * r - 1) * (2 * r)) / a1
        a4 = (2 * (r - 1) * (r - 1) * (2 * r)) / a1
        value = (Y[0] * (a3 * basisvalues[r - 1])) - a4 * basisvalues[r - 2]
        code += [L.Assign(basisvalues[r], value)]

    # Scale values
    p = L.Symbol("p")
    code += [
        L.ForRange(
            p,
            0,
            embedded_degree + 1,
            index_type=index_type,
            body=L.AssignMul(basisvalues[p], L.Sqrt(0.5 + p)))
    ]
    return code


def _generate_compute_quad_basisvalues(L, basisvalues, Y, embedded_degree,
                                       num_members):

    # Create zero-initialized array for with basisvalues
    code = [L.ArrayDecl("double", basisvalues, (num_members, ), values=0)]

    p = embedded_degree + 1
    assert p * p == num_members

    bx = [1.0]
    by = [1.0]

    if embedded_degree > 0:
        bx += [Y[0]]
        by += [Y[1]]

    # Only active if embedded_degree > 1.
    for r in range(2, p):
        a3 = (2 * r - 1) / r
        a4 = (r - 1) / r
        bx += [a3 * Y[0] * bx[r - 1] - a4 * bx[r - 2]]
        by += [a3 * Y[1] * by[r - 1] - a4 * by[r - 2]]

    for r in range(p):
        bx[r] *= numpy.sqrt(r + 0.5)
        by[r] *= numpy.sqrt(r + 0.5)

    for r in range(p * p):
        code += [L.Assign(basisvalues[r], bx[r // p] * by[r % p])]

    return code


def _generate_compute_hex_basisvalues(L, basisvalues, Y, embedded_degree,
                                      num_members):

    # Create zero-initialized array for with basisvalues
    code = [L.ArrayDecl("double", basisvalues, (num_members, ), values=0)]

    p = embedded_degree + 1
    assert p * p * p == num_members

    bx = [1.0]
    by = [1.0]
    bz = [1.0]

    if embedded_degree > 0:
        bx += [Y[0]]
        by += [Y[1]]
        bz += [Y[2]]

    # Only active if embedded_degree > 1.
    for r in range(2, p):
        a3 = (2 * r - 1) / r
        a4 = (r - 1) / r
        bx += [a3 * Y[0] * bx[r - 1] - a4 * bx[r - 2]]
        by += [a3 * Y[1] * by[r - 1] - a4 * by[r - 2]]
        bz += [a3 * Y[2] * bz[r - 1] - a4 * bz[r - 2]]

    for r in range(p):
        bx[r] *= numpy.sqrt(r + 0.5)
        by[r] *= numpy.sqrt(r + 0.5)
        bz[r] *= numpy.sqrt(r + 0.5)

    for r in range(p * p * p):
        code += [L.Assign(basisvalues[r], bx[r // (p * p)] * by[(r // p) % p] * bz[r % p])]

    return code


def _generate_compute_triangle_basisvalues(L, basisvalues, Y, embedded_degree,
                                           num_members):
    # FIAT_NEW.expansions.TriangleExpansionSet.

    def _idx2d(p, q):
        return (p + q) * (p + q + 1) // 2 + q

    # Create zero-initialized array for with basisvalues
    code = [L.ArrayDecl("double", basisvalues, (num_members, ), values=0)]

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
    code += [
        L.VariableDecl("const double", f1, (1.0 + 2.0 * Y[0] + Y[1]) / 2.0)
    ]
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
        code += [L.VariableDecl("const double", f3, f2 * f2)]
    for r in range(1, embedded_degree):
        rr = _idx2d((r + 1), 0)
        ss = _idx2d(r, 0)
        tt = _idx2d((r - 1), 0)
        A = (2 * r + 1.0) / (r + 1)
        B = r / (1.0 + r)
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
        A = 0.5 * (1 + 2 * r)
        B = 0.5 * (3 + 2 * r)
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
            A, B, C = _jrc(2 * r + 1, 0, s)
            value = (B + A * Y[1]) * basisvalues[ss] - C * basisvalues[tt]
            code += [L.Assign(basisvalues[rr], value)]

    # FIAT_NEW code (loop 4 in FIAT).
    # for p in range(n+1):
    #    for q in range(n-p+1):
    #        results[idx(p,q),:] *= math.sqrt((p+0.5)*(p+q+1.0))
    for r in range(0, embedded_degree + 1):
        for s in range(0, embedded_degree + 1 - r):
            rr = _idx2d(r, s)
            A = (r + 0.5) * (r + s + 1)
            code += [L.AssignMul(basisvalues[rr], L.Sqrt(A))]

    return code


def _generate_compute_tetrahedron_basisvalues(L, basisvalues, Y,
                                              embedded_degree, num_members):
    # FIAT_NEW.expansions.TetrahedronExpansionSet.

    # FIAT_NEW code (compute index function) TetrahedronExpansionSet.
    # def idx(p,q,r):
    #     return (p+q+r)*(p+q+r+1)*(p+q+r+2)//6 + (q+r)*(q+r+1)//2 + r
    def _idx3d(p, q, r):
        return (p + q + r) * (p + q + r + 1) * (p + q + r + 2) // 6 + (
            q + r) * (q + r + 1) // 2 + r

    # Compute helper factors.
    # FIAT_NEW code
    # factor1 = 0.5 * ( 2.0 + 2.0*x + y + z )
    # factor2 = (0.5*(y+z))**2
    # factor3 = 0.5 * ( 1 + 2.0 * y + z )
    # factor4 = 0.5 * ( 1 - z )
    # factor5 = factor4 ** 2

    # Create zero-initialized array for with basisvalues
    code = [L.ArrayDecl("double", basisvalues, (num_members, ), values=0)]

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
    code += [
        L.VariableDecl("const double", f1,
                       0.5 * (2.0 + 2.0 * Y[0] + Y[1] + Y[2]))
    ]
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
        code += [
            L.VariableDecl("const double", f2,
                           0.25 * (Y[1] + Y[2]) * (Y[1] + Y[2]))
        ]
    for r in range(1, embedded_degree):
        rr = _idx3d((r + 1), 0, 0)
        ss = _idx3d(r, 0, 0)
        tt = _idx3d((r - 1), 0, 0)
        A = (2 * r + 1.0) / (r + 1)
        B = r / (r + 1.0)
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
        term0 = 0.5 * (2.0 + 3.0 * Y[1] + Y[2])
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
        code += [
            L.VariableDecl("const double", f3, 0.5 * (1.0 + 2.0 * Y[1] + Y[2]))
        ]
        code += [L.VariableDecl("const double", f4, 0.5 * (1.0 - Y[2]))]
        code += [L.VariableDecl("const double", f5, f4 * f4)]
    for r in range(0, embedded_degree - 1):
        for s in range(1, embedded_degree - r):
            rr = _idx3d(r, (s + 1), 0)
            ss = _idx3d(r, s, 0)
            tt = _idx3d(r, s - 1, 0)
            (A, B, C) = _jrc(2 * r + 1, 0, s)
            term0 = ((A * f3) + (B * f4)) * basisvalues[ss]
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
                rr = _idx3d(r, s, (t + 1))
                ss = _idx3d(r, s, t)
                tt = _idx3d(r, s, t - 1)
                (A, B, C) = _jrc(2 * r + 2 * s + 2, 0, t)
                az_b = B + A * Y[2]
                value = (az_b * basisvalues[ss]) - (C * basisvalues[tt])
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
                A = (r + 0.5) * (r + s + 1) * (r + s + t + 1.5)
                code += [L.AssignMul(basisvalues[rr], L.Sqrt(A))]

    return code
