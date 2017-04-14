# -*- coding: utf-8 -*-
# Copyright (C) 2009-2017 Anders Logg and Martin Sandve Alnæs
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


# Note: Much of the code in this file is a direct translation
# from the old implementation in FFC, although some improvements
# have been made to the generated code.


from collections import defaultdict
import numpy
from six import string_types

from ufl import product
from ffc.uflacs.backends.ufc.generator import ufc_generator
from ffc.uflacs.backends.ufc.utils import generate_return_new_switch, generate_return_int_switch, generate_error

from ffc.uflacs.elementtables import clamp_table_small_numbers
from ffc.uflacs.backends.ufc.evaluatebasis import generate_evaluate_reference_basis
from ffc.uflacs.backends.ufc.evaluatebasis import _generate_compute_basisvalues
from ffc.uflacs.backends.ufc.evalderivs import generate_evaluate_reference_basis_derivatives
from ffc.uflacs.backends.ufc.evalderivs import _generate_combinations

# FIXME: Stop depending on legacy code
from ffc.cpp import indent
from ffc.evaluatebasis import _evaluate_basis
# from ffc.evaluatebasis import _evaluate_basis_all
from ffc.evaluatebasisderivatives import _evaluate_basis_derivatives
from ffc.evaluatebasisderivatives import _evaluate_basis_derivatives_all
# from ffc.interpolatevertexvalues import interpolate_vertex_values
from ffc.evaluatedof import evaluate_dof_and_dofs
# from ffc.evaluatedof import affine_weights

index_type = "std::size_t"

def affine_weights(dim):
    "Compute coefficents for mapping from reference to physical element"

    if dim == 1:
        return lambda x: (1.0 - x[0], x[0])
    elif dim == 2:
        return lambda x: (1.0 - x[0] - x[1], x[0], x[1])
    elif dim == 3:
        return lambda x: (1.0 - x[0] - x[1] - x[2], x[0], x[1], x[2])

def _change_variables(mapping, gdim, tdim, offset):
    """Generate code for mapping function values according to
    'mapping' and offset.

    The basics of how to map a field from a physical to the reference
    domain. (For the inverse approach -- see interpolatevertexvalues)

    Let g be a field defined on a physical domain T with physical
    coordinates x. Let T_0 be a reference domain with coordinates
    X. Assume that F: T_0 -> T such that

      x = F(X)

    Let J be the Jacobian of F, i.e J = dx/dX and let K denote the
    inverse of the Jacobian K = J^{-1}. Then we (currently) have the
    following four types of mappings:

    'affine' mapping for g:

      G(X) = g(x)

    For vector fields g:

    'contravariant piola' mapping for g:

      G(X) = det(J) K g(x)   i.e  G_i(X) = det(J) K_ij g_j(x)

    'covariant piola' mapping for g:

      G(X) = J^T g(x)          i.e  G_i(X) = J^T_ij g(x) = J_ji g_j(x)

    'double covariant piola' mapping for g:

      G(X) = J^T g(x) J     i.e. G_il(X) = J_ji g_jk(x) J_kl

    'double contravariant piola' mapping for g:

      G(X) = det(J)^2 K g(x) K^T  i.e. G_il(X)=(detJ)^2 K_ij g_jk K_lk

    """

    # meg: Various mappings must be handled both here and in
    # interpolate_vertex_values. Could this be abstracted out?

    values = L.Symbol("values")
    offset = L.Symbol("offset")

    if mapping == "affine":
        return [values[offset]]
    elif mapping == "contravariant piola":
        # Map each component from physical to reference using inverse
        # contravariant piola
        detJ = L.Symbol("detJ")
        K = L.Symbol("K")
        w = []
        for i in range(tdim):
            inner = 0.0
            for j in range(gdim):
                inner += values[j + offset]*K[i*gdim + j]
            w.append(inner*detJ)
        return w

    elif mapping == "covariant piola":
        # Map each component from physical to reference using inverse
        # covariant piola
        detJ = L.Symbol("detJ")
        J = L.Symbol("J")
        w = []
        for i in range(tdim):
            inner = 0.0
            for j in range(gdim):
                inner += values[j + offset]*J[j*tdim + i]
            w.append(inner)
        return w

    elif mapping == "double covariant piola":
        # physical to reference pullback as a covariant 2-tensor
        w = []
        J = L.Symbol("J")
        for i in range(tdim):
            for l in range(tdim):
                inner = 0.0
                for k in range(gdim):
                    for j in range(gdim):
                        inner += J[j*tdim + i] * values[j * tdim + k + offset] * J[k*tdim + l]
                w.append(inner)
        return w

    elif mapping == "double contravariant piola":
        # physical to reference using double contravariant piola
        w = []
        K = Symbol("K")
        detJ = Symbol("detJ")
        for i in range(tdim):
            for l in range(tdim):
                inner = 0.0
                for k in range(gdim):
                    for j in range(gdim):
                        inner += K[i*tdim + j] * values[j*tdim + k + offset] * K[l*tdim + k]
                w.append(inner*detJ*detJ)
        return w

    else:
        raise Exception("The mapping (%s) is not allowed" % mapping)


def jacobian(L, gdim, tdim, element_cellname):
    J = L.Symbol("J")
    coordinate_dofs = L.Symbol("coordinate_dofs")
    code = [L.Comment("Compute Jacobian"),
            L.ArrayDecl("double", J, (gdim*tdim,)),
            L.Call("compute_jacobian_"+element_cellname+"_"+str(gdim)+"d",(J, coordinate_dofs))]
    return code

def inverse_jacobian(L, gdim, tdim, element_cellname):
    K = L.Symbol("K")
    J = L.Symbol("J")
    detJ = L.Symbol("detJ")
    code = [L.Comment("Compute Inverse Jacobian and determinant"),
            L.ArrayDecl("double", K, (gdim*tdim,)),
            L.VariableDecl("double", detJ),
            L.Call("compute_jacobian_inverse_"+element_cellname+"_"+str(gdim)+"d",(K, detJ, J))]
    return code

def orientation(L):
    detJ = L.Symbol("detJ")
    cell_orientation = L.Symbol("cell_orientation")
    code = [L.Comment("Check orientation"),
            L.If(L.EQ(cell_orientation, -1),
                 [L.Throw("std::runtime_error", "cell orientation must be defined (not -1)")]),
            L.Comment("(If cell_orientation == 1 = down, multiply det(J) by -1)"),
            L.ElseIf(L.EQ(cell_orientation, 1),
                     [L.AssignMul(detJ , -1)])]
    return code

def fiat_coordinate_mapping(L, cellname, gdim):

    # Code used in evaluatebasis[|derivatives]
    x = L.Symbol("x")
    Y = L.Symbol("Y")
    coordinate_dofs = L.Symbol("coordinate_dofs")

    if cellname == "interval":
        J = L.Symbol("J")
        detJ = L.Symbol("detJ")
        if gdim == 1:
            code = [L.Comment("Get coordinates and map to the reference (FIAT) element"),
                    L.ArrayDecl("double", Y, 1, [(2*x[0] - coordinate_dofs[0] - coordinate_dofs[1])/J[0]])]
        elif gdim == 2:
            code = [L.Comment("Get coordinates and map to the reference (FIAT) element"),
                    L.ArrayDecl("double", Y, 1, [2*(L.Sqrt(L.Call("std::pow",x[0] - coordinate_dofs[0])) + L.Sqrt(L.Call("std::pow", x[1] - coordinate_dofs[1])))/detJ - 1.0])]
        elif gdim == 3:
            code = [L.Comment("Get coordinates and map to the reference (FIAT) element"),
                    L.ArrayDecl("double", Y, 1, [2*(L.Sqrt(L.Call("std::pow", (x[0] - coordinate_dofs[0], 2)) + L.Call("std::pow", (x[1] - coordinate_dofs[1], 2)) + L.Call("std::pow", (x[2] - coordinate_dofs[2], 2)))/ detJ) - 1.0])]
        else:
            error("Cannot compute interval with gdim: %d" % gdim)
    elif cellname == "triangle":
        if gdim == 2:
            C0 = L.Symbol("C0")
            C1 = L.Symbol("C1")
            J = L.Symbol("J")
            detJ = L.Symbol("detJ")
            code = [L.Comment("Compute constants"),
                    L.VariableDecl("const double", C0, coordinate_dofs[2] + coordinate_dofs[4]),
                    L.VariableDecl("const double", C1, coordinate_dofs[3] + coordinate_dofs[5]),
                    L.Comment("Get coordinates and map to the reference (FIAT) element"),
                    L.ArrayDecl("double", Y, 2, [(J[1]*(C1 - 2.0*x[1]) + J[3]*(2.0*x[0] - C0)) / detJ,
                                                 (J[0]*(2.0*x[1] - C1) + J[2]*(C0 - 2.0*x[0])) / detJ])]
        elif gdim == 3:
            K = L.Symbol("K")
            code = [L.Comment("P_FFC = J^dag (p - b), P_FIAT = 2*P_FFC - (1, 1)"),
                    L.ArrayDecl("double", Y, 2, [2*(K[0]*(x[0] - coordinate_dofs[0])
                                                    + K[1]*(x[1] - coordinate_dofs[1])
                                                    + K[2]*(x[2] - coordinate_dofs[2])) - 1.0,
                                                 2*(K[3]*(x[0] - coordinate_dofs[0])
                                                    + K[4]*(x[1] - coordinate_dofs[1])
                                                    + K[5]*(x[2] - coordinate_dofs[2])) - 1.0])]
        else:
            error("Cannot compute triangle with gdim: %d" % gdim)
    elif cellname == 'tetrahedron' and gdim == 3:
        C0 = L.Symbol("C0")
        C1 = L.Symbol("C1")
        C2 = L.Symbol("C2")
        J = L.Symbol("J")
        detJ = L.Symbol("detJ")
        d = L.Symbol("d")

        code = [L.Comment("Compute constants"),
                L.VariableDecl("const double", C0, coordinate_dofs[9]  + coordinate_dofs[6] + coordinate_dofs[3] - coordinate_dofs[0]),
                L.VariableDecl("const double", C1, coordinate_dofs[10] + coordinate_dofs[7] + coordinate_dofs[4] - coordinate_dofs[1]),
                L.VariableDecl("const double", C2, coordinate_dofs[11] + coordinate_dofs[8] + coordinate_dofs[5] - coordinate_dofs[2]),
                L.Comment("Compute subdeterminants"),
                L.ArrayDecl("const double", d, 9, [J[4]*J[8] - J[5]*J[7],
                                                   J[5]*J[6] - J[3]*J[8],
                                                   J[3]*J[7] - J[4]*J[6],
                                                   J[2]*J[7] - J[1]*J[8],
                                                   J[0]*J[8] - J[2]*J[6],
                                                   J[1]*J[6] - J[0]*J[7],
                                                   J[1]*J[5] - J[2]*J[4],
                                                   J[2]*J[3] - J[0]*J[5],
                                                   J[0]*J[4] - J[1]*J[3]]),
                L.Comment("Get coordinates and map to the reference (FIAT) element"),
                L.ArrayDecl("double", Y, 3, [(d[0]*(2.0*x[0] - C0) + d[3]*(2.0*x[1] - C1) + d[6]*(2.0*x[2] - C2)) / detJ,
                                             (d[1]*(2.0*x[0] - C0) + d[4]*(2.0*x[1] - C1) + d[7]*(2.0*x[2] - C2)) / detJ,
                                             (d[2]*(2.0*x[0] - C0) + d[5]*(2.0*x[1] - C1) + d[8]*(2.0*x[2] - C2)) / detJ])]
    else:
        error("Cannot compute %s with gdim: %d" % (cellname, gdim))

    return code

def compute_basis_values(L, data, dof_data):
    basisvalues = L.Symbol("basisvalues")
    Y = L.Symbol("Y")
    element_cellname = data["cellname"]
    embedded_degree = dof_data["embedded_degree"]
    num_members = dof_data["num_expansion_members"]
    return _generate_compute_basisvalues(L, basisvalues, Y, element_cellname, embedded_degree, num_members)

def _x_compute_basis_values(L, data, dof_data):
    # FIXME: remove this, duplicate implementation...
    # Get embedded degree.
    embedded_degree = dof_data["embedded_degree"]

    # Create zero array for basisvalues.
    # Get number of members of the expansion set.
    num_mem = dof_data["num_expansion_members"]
    code = [L.Comment("Array of basisvalues")]
    basisvalues = L.Symbol("basisvalues")
    code += [L.ArrayDecl("double", basisvalues, num_mem, 0.0)]

    # Get the element cell name
    element_cellname = data["cellname"]

    def _jrc(a, b, n):
        an = float((2 * n + 1 + a + b) * (2 * n + 2 + a + b)) / float(2 * (n + 1) * (n + 1 + a + b))
        bn = float((a * a - b * b) * (2 * n + 1 + a + b)) / float(2 * (n + 1) * (2 * n + a + b) * (n + 1 + a + b))
        cn = float((n + a) * (n + b) * (2 * n + 2 + a + b)) / float((n + 1) * (n + 1 + a + b) * (2 * n + a + b))
        return (an, bn, cn)

    # 1D
    if (element_cellname == "interval"):
        # FIAT_NEW.expansions.LineExpansionSet.
        # FIAT_NEW code
        # psitilde_as = jacobi.eval_jacobi_batch(0,0,n,ref_pts)
        # FIAT_NEW.jacobi.eval_jacobi_batch(a,b,n,xs)
        # The initial value basisvalue 0 is always 1.0
        # FIAT_NEW code
        # for ii in range(result.shape[1]):
        #    result[0,ii] = 1.0 + xs[ii,0] - xs[ii,0]
        code += [L.Comment("Compute basisvalues")]
        code += [L.Assign(basisvalues[0], 1.0)]

        # Only continue if the embedded degree is larger than zero.
        if embedded_degree > 0:

            # FIAT_NEW.jacobi.eval_jacobi_batch(a,b,n,xs).
            # result[1,:] = 0.5 * ( a - b + ( a + b + 2.0 ) * xsnew )
            # The initial value basisvalue 1 is always x
            X = L.Symbol("X")
            code += [L.Assign(basisvalues[1], X)]

            # Only active is embedded_degree > 1.
            if embedded_degree > 1:
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
                for r in range(2, embedded_degree + 1):

                    # Define helper variables
                    a1 = 2.0 * r * r * (2.0 * r - 2.0)
                    a3 = ((2.0 * r - 2.0) * (2.0 * r - 1.0) * (2.0 * r)) / a1
                    a4 = (2.0 * (r - 1.0) * (r - 1.0) * (2.0 * r)) / a1

                    code += [L.Assign(basisvalues[r], X*basisvalues[r-1]*a3 - basisvalues[r-2]*a4)]

        # Scale values.
        # FIAT_NEW.expansions.LineExpansionSet.
        # FIAT_NEW code
        # results = numpy.zeros( ( n+1 , len(pts) ) , type( pts[0][0] ) )
        # for k in range( n + 1 ):
        #    results[k,:] = psitilde_as[k,:] * math.sqrt( k + 0.5 )

        r = L.Symbol("r")
        code += [L.ForRange(r, 0, embedded_degree + 1,
                            index_type=index_type, body=[L.AssignMul(basisvalues[r], L.Sqrt(r + 0.5))])]

    # 2D
    elif (element_cellname == "triangle"):
        # FIAT_NEW.expansions.TriangleExpansionSet.

        # Compute helper factors
        # FIAT_NEW code
        # f1 = (1.0+2*x+y)/2.0
        # f2 = (1.0 - y) / 2.0
        # f3 = f2**2
        X = L.Symbol("X")
        Y = L.Symbol("Y")
        f1 = (1 + 2*X + Y)/2
        f2 = (1 - Y)/2
        f3 = f2*f2

        code += [L.Comment("Compute basisvalues")]
        # The initial value basisvalue 0 is always 1.0.
        # FIAT_NEW code
        # for ii in range( results.shape[1] ):
        #    results[0,ii] = 1.0 + apts[ii,0]-apts[ii,0]+apts[ii,1]-apts[ii,1]

        code += [L.Assign(basisvalues[0], 1.0)]

        def _idx2d(p, q):
            return (p + q) * (p + q + 1) // 2 + q

        # Only continue if the embedded degree is larger than zero.
        if embedded_degree > 0:
            # The initial value of basisfunction 1 is equal to f1.
            # FIAT_NEW code
            # results[idx(1,0),:] = f1
            code += [L.Assign(basisvalues[1], f1)]

            # NOTE: KBO: The order of the loops is VERY IMPORTANT!!

            # FIAT_NEW code (loop 1 in FIAT)
            # for p in range(1,n):
            #    a = (2.0*p+1)/(1.0+p)
            #    b = p / (p+1.0)
            #    results[idx(p+1,0)] = a * f1 * results[idx(p,0),:] \
            #        - p/(1.0+p) * f3 *results[idx(p-1,0),:]
            # FIXME: KBO: Is there an error in FIAT? why is b not used?

            # Only active is embedded_degree > 1.
            for r in range(1, embedded_degree):
                rr = _idx2d((r + 1), 0)
                ss = _idx2d(r, 0)
                tt = _idx2d(r - 1, 0)
                A = (2 * r + 1.0) / (r + 1)
                B = r / (1.0 + r)
                code += [L.Assign(basisvalues[rr],
                                  basisvalues[ss]*A*f1 - basisvalues[tt]*B*f3)]

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
                C = A + B*Y
                code += [L.Assign(basisvalues[rr],basisvalues[ss]*C)]


            # FIAT_NEW code (loop 3 in FIAT).
            # for p in range(n-1):
            #    for q in range(1,n-p):
            #        (a1,a2,a3) = jrc(2*p+1,0,q)
            #        results[idx(p,q+1),:] \
            #            = ( a1 * y + a2 ) * results[idx(p,q)] \
            #            - a3 * results[idx(p,q-1)]
            # Only active is embedded_degree > 1.
            for r in range(0, embedded_degree - 1):
                for s in range(1, embedded_degree - r):
                    rr = _idx2d(r, (s + 1))
                    ss = _idx2d(r, s)
                    tt = _idx2d(r, s - 1)
                    A, B, C = _jrc(2 * r + 1, 0, s)
                    code += [L.Assign(basisvalues[rr],
                                      basisvalues[ss]*(B + A*Y) - basisvalues[tt]*C)]

            # FIAT_NEW code (loop 4 in FIAT).
            # for p in range(n+1):
            #    for q in range(n-p+1):
            #        results[idx(p,q),:] *= math.sqrt((p+0.5)*(p+q+1.0))

            for r in range(0, embedded_degree + 1):
                for s in range(0, embedded_degree + 1 - r):
                    rr = _idx2d(r, s)
                    A = (r + 0.5) * (r + s + 1)
                    code += [L.AssignMul(basisvalues[rr], L.Sqrt(A))]

    # 3D
    elif (element_cellname == "tetrahedron"):
        # FIAT_NEW code (compute index function) TetrahedronExpansionSet.
        def _idx3d(p, q, r):
            return (p + q + r) * (p + q + r + 1) * (p + q + r + 2) // 6 + (q + r) * (q + r + 1) // 2 + r

        code += [L.Comment("Compute basisvalues")]

        # The initial value basisvalue 0 is always 1.0.
        # FIAT_NEW code
        # for ii in range( results.shape[1] ):
        #    results[0,ii] = 1.0 + apts[ii,0]-apts[ii,0]+apts[ii,1]-apts[ii,1]
        code += [L.Assign(basisvalues[0], 1.0)]

        # Only continue if the embedded degree is larger than zero.
        if embedded_degree > 0:
            X = L.Symbol("X")
            Y = L.Symbol("Y")
            Z = L.Symbol("Z")
            f1 = 0.5*(2.0 + 2.0*X + Y + Z )
            f2 = 0.5*(Y + Z)
            f2 = f2*f2
            f3 = 0.5*( 1 + 2.0 * Y + Z )
            f4 = 0.5*( 1 - Z )
            f5 = f4*f4
            # The initial value of basisfunction 1 is equal to f1.
            # FIAT_NEW code
            # results[idx(1,0),:] = f1
            code += [L.Assign(basisvalues[1], f1)]

            # NOTE: KBO: The order of the loops is VERY IMPORTANT!!
            # FIAT_NEW code (loop 1 in FIAT).
            # for p in range(1,n):
            #    a1 = ( 2.0 * p + 1.0 ) / ( p + 1.0 )
            #    a2 = p / (p + 1.0)
            #    results[idx(p+1,0,0)] = a1 * factor1 * results[idx(p,0,0)] \
            #        -a2 * factor2 * results[ idx(p-1,0,0) ]
            for r in range(1, embedded_degree):
                rr = _idx3d((r + 1), 0, 0)
                ss = _idx3d(r, 0, 0)
                tt = _idx3d((r - 1), 0, 0)
                A = (2 * r + 1.0) / (r + 1)
                B = r / (r + 1.0)
                code += [L.Assign(basisvalues[rr], A*f1*basisvalues[ss] - B*f2*basisvalues[tt])]

            # FIAT_NEW code (loop 2 in FIAT).
            # q = 1
            # for p in range(0,n):
            #    results[idx(p,1,0)] = results[idx(p,0,0)] \
            #        * ( p * (1.0 + y) + ( 2.0 + 3.0 * y + z ) / 2 )

            for r in range(0, embedded_degree):
                rr = _idx3d(r, 1, 0)
                ss = _idx3d(r, 0, 0)
                code += [L.Assign(basisvalues[rr],
                                  basisvalues[ss]*(0.5*(2 + 3*Y + Z) + r*(1 + Y)))]

            # FIAT_NEW code (loop 3 in FIAT).
            # for p in range(0,n-1):
            #    for q in range(1,n-p):
            #        (aq,bq,cq) = jrc(2*p+1,0,q)
            #        qmcoeff = aq * factor3 + bq * factor4
            #        qm1coeff = cq * factor5
            #        results[idx(p,q+1,0)] = qmcoeff * results[idx(p,q,0)] \
            #            - qm1coeff * results[idx(p,q-1,0)]

            for r in range(0, embedded_degree - 1):
                for s in range(1, embedded_degree - r):
                    rr = _idx3d(r, (s + 1), 0)
                    ss = _idx3d(r, s, 0)
                    tt = _idx3d(r, s - 1, 0)
                    (A, B, C) = _jrc(2 * r + 1, 0, s)
                    code += [L.Assign(basisvalues[rr],
                                      (A*f3 + B*f4)*basisvalues[ss] - C*f5*basisvalues[tt])]

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
                    code += [L.Assign(basisvalues[rr],
                                      basisvalues[ss]*(1 + r + s + (2 + r + s)*Z))]

            # FIAT_NEW code (loop 5 in FIAT).
            # general r by recurrence
            # for p in range(n-1):
            #     for q in range(0,n-p-1):
            #         for r in range(1,n-p-q):
            #             ar,br,cr = jrc(2*p+2*q+2,0,r)
            #             results[idx(p,q,r+1)] = \
            #                         (ar * z + br) * results[idx(p,q,r) ] \
            #                         - cr * results[idx(p,q,r-1) ]
            for r in range(embedded_degree - 1):
                for s in range(0, embedded_degree - r - 1):
                    for t in range(1, embedded_degree - r - s):
                        rr = _idx3d(r, s, (t + 1))
                        ss = _idx3d(r, s, t)
                        tt = _idx3d(r, s, t - 1)

                        (A, B, C) = _jrc(2 * r + 2 * s + 2, 0, t)
                        code += [L.Assign(basisvalues[rr],
                                          basisvalues[ss]*(A*Z + B) - basisvalues[tt]*C)]

            # FIAT_NEW code (loop 6 in FIAT).
            # for p in range(n+1):
            #    for q in range(n-p+1):
            #        for r in range(n-p-q+1):
            #            results[idx(p,q,r)] *= math.sqrt((p+0.5)*(p+q+1.0)*(p+q+r+1.5))
            for r in range(embedded_degree + 1):
                for s in range(embedded_degree - r + 1):
                    for t in range(embedded_degree - r - s + 1):
                        rr = _idx3d(r, s, t)
                        A = (r + 0.5) * (r + s + 1) * (r + s + t + 1.5)
                        code += [L.AssignMul(basisvalues[rr], L.Sqrt(A))]

    else:
        error("Cannot compute basis values for shape: %d" % element_cellname)

    return code

def tabulate_coefficients(L, dof_data):
    """This function tabulates the element coefficients that are
    generated by FIAT at compile time."""

    # Get coefficients from basis functions, computed by FIAT at compile time.
    coefficients = dof_data["coeffs"]

    # Initialise return code.
    code = [L.Comment("Table(s) of coefficients")]

    # Get number of members of the expansion set.
    num_mem = dof_data["num_expansion_members"]

    # Generate tables for each component.
    for i, coeffs in enumerate(coefficients):

        # Variable name for coefficients.
        name = L.Symbol("coefficients%d" % i)

        # Generate array of values.
        code += [L.ArrayDecl("static const double", name, num_mem, coeffs)]

    return code

def compute_values(L, data, dof_data):
    """This function computes the value of the basisfunction as the dot product
    of the coefficients and basisvalues."""

    # Initialise return code.
    code = [L.Comment("Compute value(s)")]

    # Get dof data.
    num_components = dof_data["num_components"]
    reference_offset = dof_data["reference_offset"]
    physical_offset = dof_data["physical_offset"]
    offset = reference_offset  # physical_offset # FIXME: Should be physical offset but that breaks tests

    basisvalues = L.Symbol("basisvalues")
    values = L.Symbol("values")
    r = L.Symbol("r")
    lines = []
    if data["reference_value_size"] != 1:
        # Loop number of components.
        for i in range(num_components):
            coefficients = L.Symbol("coefficients%d" % i)
            lines += [L.AssignAdd(values[i+offset], coefficients[r]*basisvalues[r])]
    else:
        coefficients = L.Symbol("coefficients0")
        lines = [L.AssignAdd(L.Dereference(values), coefficients[r]*basisvalues[r])]

    # Get number of members of the expansion set and generate loop.
    num_mem = dof_data["num_expansion_members"]
    code += [L.ForRange(r, 0, num_mem, index_type=index_type, body=lines)]

    tdim = data["topological_dimension"]
    gdim = data["geometric_dimension"]

    # Apply transformation if applicable.
    mapping = dof_data["mapping"]
    if mapping == "affine":
        pass
    elif mapping == "contravariant piola":
        code += [L.Comment("Using contravariant Piola transform to map values back to the physical element")]

        # Get temporary values before mapping.
        tmp_ref = []
        for i in range(num_components):
            tmp_ref.append(L.Symbol("tmp_ref%d" % i))
        code += [L.VariableDecl("const double", tmp_ref[i], values[i + offset])
                 for i in range(num_components)]


        # Create names for inner product.
        basis_col = [tmp_ref[j] for j in range(tdim)]
        J = L.Symbol("J")
        J = L.FlattenedArray(J, dims=(gdim, tdim))
        detJ = L.Symbol("detJ")
        for i in range(gdim):
            # Create Jacobian.
            jacobian_row = [ J[i, j] for j in range(tdim) ]

            # Create inner product and multiply by inverse of Jacobian.
            inner = 0.0
            for a,b in zip(jacobian_row, basis_col):
                inner += a*b
            value = inner/detJ
            code += [L.Assign(values[i + offset], inner/detJ)]

    elif mapping == "covariant piola":
        code += [L.Comment("Using covariant Piola transform to map values back to the physical element")]
        # Get temporary values before mapping.
        tmp_ref = []
        for i in range(num_components):
            tmp_ref.append(L.Symbol("tmp_ref%d" % i))
        code += [L.VariableDecl("const double", tmp_ref[i], values[i + offset])
                 for i in range(num_components)]

        basis_col = [tmp_ref[j] for j in range(tdim)]
        K = L.Symbol("K")
        K = L.FlattenedArray(K, dims=(tdim, gdim))
        for i in range(gdim):
            # Create inverse of Jacobian.
            inv_jacobian_column = [K[j, i] for j in range(tdim)]

            # Create inner product of basis values and inverse of Jacobian.
            inner = 0.0
            for a, b in zip(inv_jacobian_column, basis_col):
                inner += a*b
            code += [L.Assign(values[i + offset], inner)]

    elif mapping == "double covariant piola":
        code += [L.Comment("Using double covariant Piola transform to map values back to the physical element")]
        # Get temporary values before mapping.
        basis_col = []
        for i in range(num_components):
            basis_col.append(L.Symbol("tmp_ref%d" % i))
        code += [L.VariableDecl("const double", basis_col[i], values[i + offset])
                 for i in range(num_components)]

        # value = f_group(f_inner(
        #     [f_inner([f_trans("JINV", j, i, tdim, gdim, None)
        #               for j in range(tdim)],
        #              [basis_col[j * tdim + k] for j in range(tdim)])
        #      for k in range(tdim)],
        #     [f_trans("JINV", k, l, tdim, gdim, None)
        #      for k in range(tdim)]))

        K = L.Symbol("K")
        for p in range(num_components):
            # unflatten the indices
            i = p // tdim
            l = p % tdim
            # g_il = K_ji G_jk K_kl
            acc_list = []
            for k in range(tdim):
                acc = 0.0
                for j in range(tdim):
                    acc += K[j*gdim +i]*basis_col[j*tdim + k]
                acc_list.append(acc)
            inner = 0.0
            for k in range(tdim):
                inner += acc_list[k]*K[k*gdim + l]

            code += [L.Assign(values[p + offset], inner)]

    elif mapping == "double contravariant piola":
        code += [L.Comment("Using double contravariant Piola transform to map values back to the physical element")]

        # Get temporary values before mapping.
        basis_col = []
        for i in range(num_components):
            basis_col.append(L.Symbol("tmp_ref%d" % i))
        code += [L.VariableDecl("const double", basis_col[i], values[i + offset])
                 for i in range(num_components)]

        J = L.Symbol("J")
        detJ = L.Symbol("detJ")
        for p in range(num_components):
            # unflatten the indices
            i = p // tdim
            l = p % tdim

            # g_il = (det J)^(-2) Jij G_jk Jlk
            #         value = f_group(f_inner(
            #             [f_inner([f_trans("J", i, j, tdim, gdim, None)
            #                       for j in range(tdim)],
            #                      [basis_col[j * tdim + k] for j in range(tdim)])
            #              for k in range(tdim)],
            #             [f_trans("J", l, k, tdim, gdim, None) for k in range(tdim)]))

            acc_list = []
            for k in range(tdim):
                acc = 0.0
                for j in range(tdim):
                    acc += J[i*gdim + j]*basis_col[j*tdim + k]
                acc_list.append(acc)
            inner = 0.0
            for k in range(tdim):
                inner += acc_list[k]*J[l*gdim + k]

            code += [L.Assign(values[p + offset], inner/(detJ*detJ))]
    else:
        error("Unknown mapping: %s" % mapping)

    return code

def generate_element_mapping(mapping, i, num_reference_components, tdim, gdim, J, detJ, K):
    # Select transformation to apply
    if mapping == "affine":
        assert num_reference_components == 1
        num_physical_components = 1
        M_scale = 1
        M_row = [1]  # M_row[0] == 1
    elif mapping == "contravariant piola":
        assert num_reference_components == tdim
        num_physical_components = gdim
        M_scale = 1.0 / detJ
        M_row = [J[i, jj] for jj in range(tdim)]
    elif mapping == "covariant piola":
        assert num_reference_components == tdim
        num_physical_components = gdim
        M_scale = 1.0 / detJ
        M_row = [K[jj, i] for jj in range(tdim)]
    elif mapping == "double covariant piola":
        assert num_reference_components == tdim**2
        num_physical_components = gdim**2
        # g_il = K_ji G_jk K_kl = K_ji K_kl G_jk
        i0 = i // tdim  # i in the line above
        i1 = i % tdim   # l ...
        M_scale = 1.0
        M_row = [K[jj,i0]*K[kk,i1] for jj in range(tdim) for kk in range(tdim)]
    elif mapping == "double contravariant piola":
        assert num_reference_components == tdim**2
        num_physical_components = gdim**2
        # g_il = (det J)^(-2) Jij G_jk Jlk = (det J)^(-2) Jij Jlk G_jk
        i0 = i // tdim  # i in the line above
        i1 = i % tdim   # l ...
        M_scale = 1.0 / (detJ*detJ)
        M_row = [J[i0,jj]*J[i1,kk] for jj in range(tdim) for kk in range(tdim)]
    else:
        error("Unknown mapping: %s" % mapping)
    return M_scale, M_row, num_physical_components


class ufc_finite_element(ufc_generator):
    "Each function maps to a keyword in the template. See documentation of ufc_generator."
    def __init__(self):
        ufc_generator.__init__(self, "finite_element")

    def cell_shape(self, L, cell_shape):
        return L.Return(L.Symbol("ufc::shape::" + cell_shape))

    def topological_dimension(self, L, topological_dimension):
        return L.Return(topological_dimension)

    def geometric_dimension(self, L, geometric_dimension):
        return L.Return(geometric_dimension)

    def space_dimension(self, L, space_dimension):
        return L.Return(space_dimension)

    def value_rank(self, L, value_shape):
        return L.Return(len(value_shape))

    def value_dimension(self, L, value_shape):
        return generate_return_int_switch(L, "i", value_shape, 1)

    def value_size(self, L, value_shape):
        return L.Return(product(value_shape))

    def reference_value_rank(self, L, reference_value_shape):
        return L.Return(len(reference_value_shape))

    def reference_value_dimension(self, L, reference_value_shape):
        return generate_return_int_switch(L, "i", reference_value_shape, 1)

    def reference_value_size(self, L, reference_value_shape):
        return L.Return(product(reference_value_shape))

    def degree(self, L, degree):
        return L.Return(degree)

    def family(self, L, family):
        return L.Return(L.LiteralString(family))

    def num_sub_elements(self, L, num_sub_elements):
        return L.Return(num_sub_elements)

    def create_sub_element(self, L, ir):
        classnames = ir["create_sub_element"]
        return generate_return_new_switch(L, "i", classnames, factory=ir["jit"])

    def evaluate_basis(self, L, ir, parameters):
        legacy_code = indent(_evaluate_basis(ir["evaluate_basis"]), 4)
        #        print(legacy_code)

        data = ir["evaluate_basis"]
        if isinstance(data, string_types):
            return format["exception"]("evaluate_basis: %s" % data)

        # Get the element cell name and geometric dimension.
        element_cellname = data["cellname"]
        gdim = data["geometric_dimension"]
        tdim = data["topological_dimension"]

        # Generate run time code to evaluate an element basisfunction at an
        # arbitrary point. The value(s) of the basisfunction is/are
        # computed as in FIAT as the dot product of the coefficients (computed at compile time)
        # and basisvalues which are dependent on the coordinate and thus have to be computed at
        # run time.

        # The function should work for all elements supported by FIAT, but it remains
        # untested for tensor valued elements.

        # Get code snippets for Jacobian, Inverse of Jacobian and mapping of
        # coordinates from physical element to the FIAT reference element.

        code = jacobian(L, gdim, tdim, element_cellname)
        code += inverse_jacobian(L, gdim, tdim, element_cellname)
        if data["needs_oriented"]:
            code += orientation(L)

        need_mp = False
        for  dof_data in data["dofs_data"]:
            if dof_data['embedded_degree'] > 0 :
                need_mp = True
        if need_mp:
            code += fiat_coordinate_mapping(L, element_cellname, gdim)

        reference_value_size = data["reference_value_size"]
        code += [L.Comment("Reset values")]
        dof_values = L.Symbol("values")
        if reference_value_size == 1:
            # Reset values as a pointer.
            code += [L.Assign(L.Dereference(dof_values), 0.0)]
        else:
            code += [L.MemZero(dof_values, reference_value_size)]

        # Create code for all basis values (dofs).
        dof_cases = []
        for f, dof_data in enumerate(data["dofs_data"]):
            dof_code = compute_basis_values(L, data, dof_data)
            dof_code += tabulate_coefficients(L, dof_data)
            dof_code += compute_values(L, data, dof_data)
            dof_cases.append((f, dof_code))

        code += [L.Switch(L.Symbol("i"), dof_cases)]

        #        print(L.StatementList(code))
        return code

    def evaluate_basis_all(self, L, ir, parameters):

        data=ir["evaluate_basis"]
        physical_value_size = data["physical_value_size"]
        space_dimension = data["space_dimension"]

        x = L.Symbol("x")
        coordinate_dofs = L.Symbol("coordinate_dofs")
        cell_orientation = L.Symbol("cell_orientation")
        values = L.Symbol("values")

        # Special case where space dimension is one (constant elements).
        if space_dimension == 1:
            code = [L.Comment("Element is constant, calling evaluate_basis."),
                    L.Call("evaluate_basis",
                           (0, values, x, coordinate_dofs, cell_orientation))]
            return code

        r = L.Symbol("r")
        dof_values = L.Symbol("dof_values")
        if physical_value_size == 1:
            code = [ L.Comment("Helper variable to hold value of a single dof."),
                     L.VariableDecl("double", dof_values, 0.0),
                     L.Comment("Loop dofs and call evaluate_basis"),
                     L.ForRange(r, 0, space_dimension, index_type=index_type,
                                body=[L.Call("evaluate_basis",
                                             (r, L.AddressOf(dof_values), x,
                                              coordinate_dofs, cell_orientation)),
                                      L.Assign(values[r], dof_values)]
                               )
                   ]
        else:
            s = L.Symbol("s")
            code = [L.Comment("Helper variable to hold values of a single dof."),
                    L.ArrayDecl("double", dof_values, physical_value_size, 0.0),
                    L.Comment("Loop dofs and call evaluate_basis"),
                    L.ForRange(r, 0, space_dimension, index_type=index_type,
                               body=[L.Call("evaluate_basis",
                                             (r, dof_values, x,
                                              coordinate_dofs, cell_orientation)),
                                     L.ForRange(s, 0, physical_value_size,
                                                index_type=index_type,
                                                body=[L.Assign(values[r*physical_value_size+s], dof_values[s])])
                                    ]
                              )
                   ]

        return code

    def evaluate_basis_derivatives(self, L, ir, parameters):
        # FIXME: Get rid of this
        # FIXME: port this
        legacy_code = indent(_evaluate_basis_derivatives(ir["evaluate_basis"]), 4)
        return legacy_code


    def evaluate_basis_derivatives_all(self, L, ir, parameters):
        # FIXME: port this
        use_legacy = 1
        if use_legacy:
            return indent(_evaluate_basis_derivatives_all(ir["evaluate_basis"]), 4)

        """
        // Legacy version:
        evaluate_basis_derivatives_all(std::size_t n,
                                       double * values,
                                       const double * x,
                                       const double * coordinate_dofs,
                                       int cell_orientation)
        // Suggestion for new version:
        new_evaluate_basis_derivatives(double * values,
                                       std::size_t order,
                                       std::size_t num_points,
                                       const double * x,
                                       const double * coordinate_dofs,
                                       int cell_orientation,
                                       const ufc::coordinate_mapping * cm)
        """

        # TODO: This is a refactoring step to allow rewriting code
        # generation to use coordinate_mapping in one stage, before
        # making it available as an argument from dolfin in the next stage.
        affine_coordinate_mapping_classname = ir["affine_coordinate_mapping_classname"]

        # Output arguments:
        values = L.Symbol("values")

        # Input arguments:
        #order = L.Symbol("order")
        order = L.Symbol("n")
        x = L.Symbol("x")
        coordinate_dofs = L.Symbol("coordinate_dofs")
        cell_orientation = L.Symbol("cell_orientation")

        # Internal variables:
        #num_points = L.Symbol("num_points")
        num_points = 1  # Always 1 in legacy API
        reference_values = L.Symbol("reference_values")
        X = L.Symbol("X")
        J = L.Symbol("J")
        detJ = L.Symbol("detJ")
        K = L.Symbol("K")
        ip = L.Symbol("ip")

        gdim = ir["geometric_dimension"]
        tdim = ir["topological_dimension"]

        code = [
            # Create local affine coordinate mapping object
            # TODO: Get this as input instead to support non-affine
            L.VariableDecl(affine_coordinate_mapping_classname, "cm"),
            L.ForRange(ip, 0, num_points, index_type=index_type, body=[
                L.ArrayDecl("double", X, (tdim,)),
                L.ArrayDecl("double", J, (gdim*tdim,)),
                L.ArrayDecl("double", detJ, (1,)),
                L.ArrayDecl("double", K, (tdim*gdim,)),
                L.Call("cm.compute_reference_geometry",
                       (X, J, detJ, K, num_points, x, coordinate_dofs, cell_orientation)),
                L.Call("evaluate_reference_basis_derivatives",
                       (reference_values, order, num_points, X)),
                L.Call("transform_reference_basis_derivatives",
                       (values, order, num_points, reference_values, X, J, detJ, K, cell_orientation)),
            ])
        ]
        return code

    def evaluate_dof(self, L, ir, parameters):
        # FIXME: Get rid of this
        # FIXME: port this
        use_legacy = 1
        if use_legacy:
            # Codes generated together
            (evaluate_dof_code, evaluate_dofs_code) \
              = evaluate_dof_and_dofs(ir["evaluate_dof"])
            return indent(evaluate_dof_code, 4)

    def evaluate_dofs(self, L, ir, parameters):
        """Generate code for evaluate_dofs."""
        """
        - evaluate_dof needs to be split into invert_mapping + evaluate_dof or similar?

          f = M fhat;  nu(f) = nu(M fhat) = nuhat(M^-1 f) = sum_i w_i M^-1 f(x_i)

          // Get fixed set of points on reference element
          num_points = element->num_dof_evaluation_points();
          double X[num_points*tdim];
          element->tabulate_dof_evaluation_points(X);

          // Compute geometry in these points
          domain->compute_geometry(reference_points, num_point, X, J, detJ, K, coordinate_dofs, cell_orientation);

          // Computed by dolfin
          for ip
            fvalues[ip][:] = f.evaluate(point[ip])[:];

          // Finally: nu_j(f) = sum_component sum_ip weights[j][ip][component] fvalues[ip][component]
          element->evaluate_dofs(fdofs, fvalues, J, detJ, K)
        """
        # FIXME: port this, then translate into reference version
        use_legacy = 1
        if use_legacy:
            # Codes generated together
            (evaluate_dof_code, evaluate_dofs_code) \
              = evaluate_dof_and_dofs(ir["evaluate_dof"])
            return indent(evaluate_dofs_code, 4)

    def interpolate_vertex_values(self, L, ir, parameters):
        # legacy_code = indent(interpolate_vertex_values(ir["interpolate_vertex_values"]), 4)
        #        print(legacy_code)

        irdata = ir["interpolate_vertex_values"]
        # Raise error if interpolate_vertex_values is ill-defined
        if not irdata:
            msg = "interpolate_vertex_values is not defined for this element"
            return generate_error(L, msg, parameters["convert_exceptions_to_warnings"])

        # Add code for Jacobian if necessary
        code = []
        gdim = irdata["geometric_dimension"]
        tdim = irdata["topological_dimension"]
        element_cellname = ir["evaluate_basis"]["cellname"]
        if irdata["needs_jacobian"]:
            code += jacobian(L, gdim, tdim, element_cellname)
            code += inverse_jacobian(L, gdim, tdim, element_cellname)
            if irdata["needs_oriented"] and tdim != gdim:
                code += orientation(L)

        # Compute total value dimension for (mixed) element
        total_dim = irdata["physical_value_size"]

        # Generate code for each element
        value_offset = 0
        space_offset = 0
        for data in irdata["element_data"]:
            # Add vertex interpolation for this element
            code += [L.Comment("Evaluate function and change variables")]

            # Extract vertex values for all basis functions
            vertex_values = data["basis_values"]
            value_size = data["physical_value_size"]
            space_dim = data["space_dim"]
            mapping = data["mapping"]

            # Create code for each value dimension:
            for k in range(value_size):
                # Create code for each vertex x_j
                for (j, values_at_vertex) in enumerate(vertex_values):

                    if value_size == 1:
                        values_at_vertex = [values_at_vertex]

                    values = clamp_table_small_numbers(values_at_vertex)

                    # Map basis functions using appropriate mapping
                    # FIXME: sort out all non-affine mappings and make into a function
                    # components = change_of_variables(values_at_vertex, k)

                    if mapping == 'affine':
                        w = values[k]
                    elif mapping == 'contravariant piola':
                        detJ = L.Symbol("detJ")
                        J = L.Symbol("J")
                        w = []
                        for index in range(space_dim):
                            inner = 0.0
                            for p in range(tdim):
                                inner += J[p+k*tdim]*values[p][index]
                            w.append(inner/detJ)
                    elif mapping == 'covariant piola':
                        K = L.Symbol("K")
                        w = []
                        for index in range(space_dim):
                            acc_sum = 0.0
                            for p in range(tdim):
                                acc_sum += K[k+p*gdim]*values[p][index]
                            w.append(acc_sum)
                    elif mapping == 'double covariant piola':
                        K = L.Symbol("K")
                        w = []
                        for index in range(space_dim):
                            acc_sum = 0.0
                            for p in range(tdim):
                                for q in range(tdim):
                                    acc_sum += K[k//tdim + p*gdim]*values[p][q][index]*K[k % tdim + q*gdim]
                            w.append(acc_sum)
                    elif mapping == 'double contravariant piola':
                        J = L.Symbol("J")
                        detJ = L.Symbol("detJ")
                        w = []
                        for index in range(space_dim):
                            acc_sum = 0.0
                            for p in range(tdim):
                                for q in range(tdim):
                                    acc_sum += J[p + (k//tdim)*tdim]*values[p][q][index]*J[q + (k % tdim)*tdim]
                            acc_sum /= (detJ*detJ)
                            w.append(acc_sum)
                    else:
                        raise RuntimeError("Mapping not implemented")

                    # Contract coefficients and basis functions
                    dof_values = L.Symbol("dof_values")
                    dof_list = [dof_values[i + space_offset] for i in range(space_dim)]
                    acc_value = 0.0
                    for p, q in zip(dof_list, w):
                        acc_value += p*q

                    # Assign value to correct vertex
                    index = j * total_dim + (k + value_offset)
                    v_values = L.Symbol("vertex_values")
                    code += [L.Assign(v_values[index], acc_value)]

            # Update offsets for value- and space dimension
            value_offset += data["physical_value_size"]
            space_offset += data["space_dim"]

            #        print(L.StatementList(code))
        return code

    def tabulate_dof_coordinates(self, L, ir, parameters):
        ir = ir["tabulate_dof_coordinates"]

        # Raise error if tabulate_dof_coordinates is ill-defined
        if not ir:
            msg = "tabulate_dof_coordinates is not defined for this element"
            return generate_error(L, msg, parameters["convert_exceptions_to_warnings"])

        # Extract coordinates and cell dimension
        gdim = ir["gdim"]
        tdim = ir["tdim"]
        points = ir["points"]

        # Output argument
        dof_coordinates = L.FlattenedArray(L.Symbol("dof_coordinates"),
                                           dims=(len(points), gdim))

        # Input argument
        coordinate_dofs = L.Symbol("coordinate_dofs")

        # Loop indices
        i = L.Symbol("i")
        k = L.Symbol("k")
        ip = L.Symbol("ip")

        # Basis symbol
        phi = L.Symbol("phi")

        # TODO: Get rid of all places that use affine_weights, assumes affine mesh
        # Create code for evaluating affine coordinate basis functions
        num_scalar_xdofs = tdim + 1
        cg1_basis = affine_weights(tdim)
        phi_values = numpy.asarray([phi_comp for X in points for phi_comp in cg1_basis(X)])
        assert len(phi_values) == len(points) * num_scalar_xdofs

        # TODO: Use precision parameter here
        phi_values = clamp_table_small_numbers(phi_values)

        code = [
            L.Assign(
                dof_coordinates[ip][i],
                sum(phi_values[ip*num_scalar_xdofs + k] * coordinate_dofs[gdim*k + i]
                    for k in range(num_scalar_xdofs))
            )
            for ip in range(len(points))
            for i in range(gdim)
        ]

        # FIXME: This code assumes an affine coordinate field.
        #        To get around that limitation, make this function take another argument
        #            const ufc::coordinate_mapping * cm
        #        and generate code like this:
        """
        index_type X[tdim*num_dofs];
        tabulate_dof_coordinates(X);
        cm->compute_physical_coordinates(x, X, coordinate_dofs);
        """

        return code

    def tabulate_reference_dof_coordinates(self, L, ir, parameters):
        # TODO: Change signature to avoid copy? E.g.
        # virtual const std::vector<double> & tabulate_reference_dof_coordinates() const = 0;
        # See integral::enabled_coefficients for example

        # TODO: ensure points is a numpy array,
        #   get tdim from points.shape[1],
        #   place points in ir directly instead of the subdict
        ir = ir["tabulate_dof_coordinates"]

        # Raise error if tabulate_reference_dof_coordinates is ill-defined
        if not ir:
            msg = "tabulate_reference_dof_coordinates is not defined for this element"
            return generate_error(L, msg, parameters["convert_exceptions_to_warnings"])

        # Extract coordinates and cell dimension
        tdim = ir["tdim"]
        points = ir["points"]

        # Output argument
        reference_dof_coordinates = L.Symbol("reference_dof_coordinates")

        # Reference coordinates
        dof_X = L.Symbol("dof_X")
        dof_X_values = [X[jj] for X in points for jj in range(tdim)]
        decl = L.ArrayDecl("static const double", dof_X,
                           (len(points) * tdim,), values=dof_X_values)
        copy = L.MemCopy(dof_X, reference_dof_coordinates, tdim*len(points))

        code = [decl, copy]
        return code

    def evaluate_reference_basis(self, L, ir, parameters):
        data = ir["evaluate_basis"]
        if isinstance(data, string_types):
            msg = "evaluate_reference_basis: %s" % data
            return generate_error(L, msg, parameters["convert_exceptions_to_warnings"])

        return generate_evaluate_reference_basis(L, data, parameters)

    def evaluate_reference_basis_derivatives(self, L, ir, parameters):
        data = ir["evaluate_basis"]
        if isinstance(data, string_types):
            msg = "evaluate_reference_basis_derivatives: %s" % data
            return generate_error(L, msg, parameters["convert_exceptions_to_warnings"])

        return generate_evaluate_reference_basis_derivatives(L, data, parameters)

    def transform_reference_basis_derivatives(self, L, ir, parameters):
        data = ir["evaluate_basis"]
        if isinstance(data, string_types):
            msg = "transform_reference_basis_derivatives: %s" % data
            return generate_error(L, msg, parameters["convert_exceptions_to_warnings"])

        # Get some known dimensions
        #element_cellname = data["cellname"]
        gdim = data["geometric_dimension"]
        tdim = data["topological_dimension"]
        max_degree = data["max_degree"]
        reference_value_size = data["reference_value_size"]
        physical_value_size = data["physical_value_size"]
        num_dofs = len(data["dofs_data"])

        max_g_d = gdim**max_degree
        max_t_d = tdim**max_degree

        # Output arguments
        values_symbol = L.Symbol("values")

        # Input arguments
        order = L.Symbol("order")
        num_points = L.Symbol("num_points")  # FIXME: Currently assuming 1 point?
        reference_values = L.Symbol("reference_values")
        J = L.Symbol("J")
        detJ = L.Symbol("detJ")
        K = L.Symbol("K")

        # Internal variables
        transform = L.Symbol("transform")

        # Indices, I've tried to use these for a consistent purpose
        ip = L.Symbol("ip") # point
        i = L.Symbol("i")   # physical component
        j = L.Symbol("j")   # reference component
        k = L.Symbol("k")   # order
        r = L.Symbol("r")   # physical derivative number
        s = L.Symbol("s")   # reference derivative number
        d = L.Symbol("d")   # dof

        combinations_code = []
        if max_degree == 0:
            # Don't need combinations
            num_derivatives_t = 1  # TODO: I think this is the right thing to do to make this still work for order=0?
            num_derivatives_g = 1
        elif tdim == gdim:
            num_derivatives_t = L.Symbol("num_derivatives")
            num_derivatives_g = num_derivatives_t
            combinations_code += [
                L.VariableDecl("const " + index_type, num_derivatives_t,
                               L.Call("std::pow", (tdim, order))),
            ]

            # Add array declarations of combinations
            combinations_code_t, combinations_t = _generate_combinations(L, tdim, max_degree, order, num_derivatives_t)
            combinations_code += combinations_code_t
            combinations_g = combinations_t
        else:
            num_derivatives_t = L.Symbol("num_derivatives_t")
            num_derivatives_g = L.Symbol("num_derivatives_g")
            combinations_code += [
                L.VariableDecl("const " + index_type, num_derivatives_t,
                               L.Call("std::pow", (tdim, order))),
                L.VariableDecl("const " + index_type, num_derivatives_g,
                               L.Call("std::pow", (gdim, order))),
            ]
            # Add array declarations of combinations
            combinations_code_t, combinations_t = _generate_combinations(L, tdim, max_degree, order, num_derivatives_t, suffix="_t")
            combinations_code_g, combinations_g = _generate_combinations(L, gdim, max_degree, order, num_derivatives_g, suffix="_g")
            combinations_code += combinations_code_t
            combinations_code += combinations_code_g

        # Define expected dimensions of argument arrays
        J = L.FlattenedArray(J, dims=(num_points, gdim, tdim))
        detJ = L.FlattenedArray(detJ, dims=(num_points,))
        K = L.FlattenedArray(K, dims=(num_points, tdim, gdim))

        values = L.FlattenedArray(values_symbol,
            dims=(num_points, num_dofs, num_derivatives_g, physical_value_size))
        reference_values = L.FlattenedArray(reference_values,
            dims=(num_points, num_dofs, num_derivatives_t, reference_value_size))

        # Generate code to compute the derivative transform matrix
        transform_matrix_code = [
            # Initialize transform matrix to all 1.0
            L.ArrayDecl("double", transform, (max_g_d, max_t_d)),
            L.ForRanges(
                (r, 0, num_derivatives_g),
                (s, 0, num_derivatives_t),
                index_type=index_type,
                body=L.Assign(transform[r, s], 1.0)
            ),
            ]
        if max_degree > 0:
            transform_matrix_code += [
                # Compute transform matrix entries, each a product of K entries
                L.ForRanges(
                    (r, 0, num_derivatives_g),
                    (s, 0, num_derivatives_t),
                    (k, 0, order),
                    index_type=index_type,
                    body=L.AssignMul(transform[r, s],
                                     K[ip, combinations_t[s, k], combinations_g[r, k]])
                ),
            ]

        # Initialize values to 0, will be added to inside loops
        values_init_code = [
            L.MemZero(values_symbol, num_points * num_dofs * num_derivatives_g * physical_value_size),
            ]

        # Make offsets available in generated code
        reference_offsets = L.Symbol("reference_offsets")
        physical_offsets = L.Symbol("physical_offsets")
        dof_attributes_code = [
            L.ArrayDecl("const " + index_type, reference_offsets, (num_dofs,),
                        values=[dof_data["reference_offset"] for dof_data in data["dofs_data"]]),
            L.ArrayDecl("const " + index_type, physical_offsets, (num_dofs,),
                        values=[dof_data["physical_offset"] for dof_data in data["dofs_data"]]),
            ]

        # Build dof lists for each mapping type
        mapping_dofs = defaultdict(list)
        for idof, dof_data in enumerate(data["dofs_data"]):
            mapping_dofs[dof_data["mapping"]].append(idof)

        # Generate code for each mapping type
        d = L.Symbol("d")
        transform_apply_code = []
        for mapping in sorted(mapping_dofs):
            # Get list of dofs using this mapping
            idofs = mapping_dofs[mapping]

            # Select iteration approach over dofs
            if idofs == list(range(idofs[0], idofs[-1]+1)):
                # Contiguous
                dofrange = (d, idofs[0], idofs[-1]+1)
                idof = d
            else:
                # Stored const array of dof indices
                idofs_symbol = L.Symbol("%s_dofs" % mapping.replace(" ", "_"))
                dof_attributes_code += [
                    L.ArrayDecl("const " + index_type, idofs_symbol,
                                (len(idofs),), values=idofs),
                ]
                dofrange = (d, 0, len(idofs))
                idof = idofs_symbol[d]

            # NB! Array access to offsets, these are not Python integers
            reference_offset = reference_offsets[idof]
            physical_offset = physical_offsets[idof]

            # How many components does each basis function with this mapping have?
            # This should be uniform, i.e. there should be only one element in this set:
            num_reference_components, = set(data["dofs_data"][i]["num_components"] for i in idofs)

            M_scale, M_row, num_physical_components = generate_element_mapping(
                mapping, i,
                num_reference_components, tdim, gdim,
                J[ip], detJ[ip], K[ip]
            )

            transform_apply_body = [
                L.AssignAdd(values[ip, idof, r, physical_offset + k],
                            transform[r, s] * reference_values[ip, idof, s, reference_offset + k])
                for k in range(num_physical_components)
            ]

            msg = "Using %s transform to map values back to the physical element." % mapping.replace("piola", "Piola")

            mapped_value = L.Symbol("mapped_value")
            transform_apply_code += [
                L.ForRanges(
                    dofrange,
                    (s, 0, num_derivatives_t),
                    (i, 0, num_physical_components),
                    index_type=index_type, body=[
                        # Unrolled application of mapping to one physical component,
                        # for affine this automatically reduces to
                        #   mapped_value = reference_values[..., reference_offset]
                        L.Comment(msg),
                        L.VariableDecl("const double", mapped_value,
                                       M_scale * sum(M_row[jj] * reference_values[ip, idof, s, reference_offset + jj]
                                                     for jj in range(num_reference_components))),
                        # Apply derivative transformation, for order=0 this reduces to
                        # values[ip,idof,0,physical_offset+i] = transform[0,0]*mapped_value
                        L.Comment("Mapping derivatives back to the physical element"),
                        L.ForRanges(
                            (r, 0, num_derivatives_g),
                            index_type=index_type, body=[
                                L.AssignAdd(values[ip, idof, r, physical_offset + i],
                                            transform[r, s] * mapped_value)
                        ])
                ])
            ]

        # Transform for each point
        point_loop_code = [
            L.ForRange(ip, 0, num_points, index_type=index_type, body=(
                transform_matrix_code
                + transform_apply_code
            ))
        ]

        # Join code
        code = (
            combinations_code
            + values_init_code
            + dof_attributes_code
            + point_loop_code
        )
        return code
