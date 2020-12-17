# Copyright (C) 2020 Matthew Scroggs
#
# This file is part of FFCX.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import cffi
import numpy as np
import pytest

import ffcx.codegeneration.jit
import ufl
import sympy


def float_to_type(name):
    """Map a string name to C and NumPy types"""
    if name == "double":
        return "double", np.float64
    elif name == "double complex":
        return "double _Complex", np.complex128
    elif name == "float":
        return "float", np.float32
    elif name == "float complex":
        return "float _Complex", np.complex64
    elif name == "long double":
        return "long double", np.longdouble
    else:
        raise RuntimeError("Unknown C type for: {}".format(name))


def lagrange_triangle_symbolic(order, corners=[(1, 0), (2, 0), (0, 1)], fun=lambda i: i):
    x = sympy.Symbol("x")
    y = sympy.Symbol("y")
    from sympy import S
    poly_basis = [x**i * y**j for i in range(order + 1) for j in range(order + 1 - i)]
    # vertices
    eval_points = [S(c) for c in corners]
    # edges
    for e in [(1, 2), (0, 2), (0, 1)]:
        p0 = corners[e[0]]
        p1 = corners[e[1]]
        eval_points += [tuple(S(a) + sympy.Rational((b - a) * i, order) for a, b in zip(p0, p1))
                        for i in range(1, order)]
    # face
    for f in [(0, 1, 2)]:
        p0 = corners[f[0]]
        p1 = corners[f[1]]
        p2 = corners[f[2]]
        eval_points += [tuple(S(a) + sympy.Rational((b - a) * i, order)
                        + sympy.Rational((c - a) * j, order) for a, b, c in zip(p0, p1, p2))
                        for i in range(1, order) for j in range(1, order - i)]

    dual_mat = [[f.subs(x, p[0]).subs(y, p[1]) for p in eval_points] for f in poly_basis]
    dual_mat = sympy.Matrix(dual_mat)
    mat = dual_mat.inv()
    functions = [sum(i * j for i, j in zip(mat.row(k), poly_basis)) for k in range(mat.rows)]
    results = []
    for f in functions:
        integrand = fun(f)
        results.append(integrand.integrate((x, 1 - y, 2 - 2 * y), (y, 0, 1)))
    return results


@pytest.mark.parametrize("mode", ["double"])
@pytest.mark.parametrize("sym_fun,ufl_fun", [
    (lambda i: i, lambda i: i),
    (lambda i: i.diff("x"), lambda i: ufl.grad(i)[0]),
    (lambda i: i.diff("y"), lambda i: ufl.grad(i)[1])])
@pytest.mark.parametrize("order", [1, 2, 3, 4, 5])
def test_lagrange_triangle(compile_args, order, mode, sym_fun, ufl_fun):
    sym = lagrange_triangle_symbolic(order, fun=sym_fun)
    cell = ufl.triangle
    element = ufl.FiniteElement("Lagrange", cell, order)
    v = ufl.TestFunction(element)

    a = ufl_fun(v) * ufl.dx
    forms = [a]
    compiled_forms, module = ffcx.codegeneration.jit.compile_forms(
        forms, parameters={'scalar_type': mode}, cffi_extra_compile_args=compile_args)

    ffi = cffi.FFI()
    form0 = compiled_forms[0][0]

    assert form0.num_cell_integrals == 1
    ids = np.zeros(form0.num_cell_integrals, dtype=np.intc)
    form0.get_cell_integral_ids(ffi.cast('int *', ids.ctypes.data))
    assert ids[0] == -1

    default_integral = form0.create_cell_integral(ids[0])

    c_type, np_type = float_to_type(mode)
    b = np.zeros((order + 2) * (order + 1) // 2, dtype=np_type)
    w = np.array([], dtype=np_type)

    coords = np.array([1.0, 0.0, 2.0, 0.0, 0.0, 1.0], dtype=np.float64)
    default_integral.tabulate_tensor(
        ffi.cast('{type} *'.format(type=c_type), b.ctypes.data),
        ffi.cast('{type} *'.format(type=c_type), w.ctypes.data),
        ffi.NULL,
        ffi.cast('double *', coords.ctypes.data), ffi.NULL, ffi.NULL, 0)

    # Check that the result is the same as for sympy
    assert np.allclose(b, [float(i) for i in sym])

    # Check that passing in permutations correctly reverses edges
    for perm in range(1, 8):
        perm_b = np.zeros((order + 2) * (order + 1) // 2, dtype=np_type)
        default_integral.tabulate_tensor(
            ffi.cast('{type} *'.format(type=c_type), perm_b.ctypes.data),
            ffi.cast('{type} *'.format(type=c_type), w.ctypes.data),
            ffi.NULL,
            ffi.cast('double *', coords.ctypes.data), ffi.NULL, ffi.NULL, perm)

        for edge in range(3):
            start = 3 + (order - 1) * edge
            end = start + order - 1
            if perm >> edge & 1 == 1:
                assert np.allclose(b[start: end], perm_b[end - 1: start - 1: -1])
            else:
                assert np.allclose(b[start: end], perm_b[start: end])


def lagrange_tetrahedron_symbolic(order, corners=[(1, 0, 0), (2, 0, 0), (0, 1, 0), (0, 0, 1)], fun=lambda i: i):
    x = sympy.Symbol("x")
    y = sympy.Symbol("y")
    z = sympy.Symbol("z")
    from sympy import S
    poly_basis = [
        x**i * y**j * z**k for i in range(order + 1) for j in range(order + 1 - i)
        for k in range(order + 1 - i - j)]
    # vertices
    eval_points = [S(c) for c in corners]
    # edges
    for e in [(2, 3), (1, 3), (1, 2), (0, 3), (0, 2), (0, 1)]:
        p0 = corners[e[0]]
        p1 = corners[e[1]]
        eval_points += [tuple(S(a) + sympy.Rational((b - a) * i, order) for a, b in zip(p0, p1))
                        for i in range(1, order)]
    # face
    for f in [(1, 2, 3), (0, 2, 3), (0, 1, 3), (0, 1, 2)]:
        p0 = corners[f[0]]
        p1 = corners[f[1]]
        p2 = corners[f[2]]
        eval_points += [tuple(S(a) + sympy.Rational((b - a) * i, order)
                        + sympy.Rational((c - a) * j, order) for a, b, c in zip(p0, p1, p2))
                        for i in range(1, order) for j in range(1, order - i)]
    # interior
    for v in [(0, 1, 2, 3)]:
        p0 = corners[v[0]]
        p1 = corners[v[1]]
        p2 = corners[v[2]]
        p3 = corners[v[3]]
        eval_points += [tuple(S(a) + sympy.Rational((b - a) * i, order) + sympy.Rational((c - a) * j, order)
                        + sympy.Rational((d - a) * k, order) for a, b, c, d in zip(p0, p1, p2, p3))
                        for i in range(1, order) for j in range(1, order - i) for k in range(1, order - i - j)]

    dual_mat = [[f.subs(x, p[0]).subs(y, p[1]).subs(z, p[2]) for p in eval_points] for f in poly_basis]
    dual_mat = sympy.Matrix(dual_mat)
    mat = dual_mat.inv()
    functions = [sum(i * j for i, j in zip(mat.row(k), poly_basis)) for k in range(mat.rows)]
    results = []
    for f in functions:
        integrand = fun(f)
        results.append(integrand.integrate((x, 1 - y - z, 2 - 2 * y - 2 * z), (y, 0, 1 - z), (z, 0, 1)))
    return results


@pytest.mark.parametrize("mode", ["double"])
@pytest.mark.parametrize("sym_fun,ufl_fun", [
    (lambda i: i, lambda i: i),
    (lambda i: i.diff("x"), lambda i: ufl.grad(i)[0]),
    (lambda i: i.diff("y"), lambda i: ufl.grad(i)[1])])
@pytest.mark.parametrize("order", [1, 2, 3, 4])
def test_lagrange_tetrahedron(compile_args, order, mode, sym_fun, ufl_fun):
    sym = lagrange_tetrahedron_symbolic(order, fun=sym_fun)
    cell = ufl.tetrahedron
    element = ufl.FiniteElement("Lagrange", cell, order)
    v = ufl.TestFunction(element)

    a = ufl_fun(v) * ufl.dx
    forms = [a]
    compiled_forms, module = ffcx.codegeneration.jit.compile_forms(
        forms, parameters={'scalar_type': mode}, cffi_extra_compile_args=compile_args)

    ffi = cffi.FFI()
    form0 = compiled_forms[0][0]

    assert form0.num_cell_integrals == 1
    ids = np.zeros(form0.num_cell_integrals, dtype=np.intc)
    form0.get_cell_integral_ids(ffi.cast('int *', ids.ctypes.data))
    assert ids[0] == -1

    default_integral = form0.create_cell_integral(ids[0])

    c_type, np_type = float_to_type(mode)
    b = np.zeros((order + 3) * (order + 2) * (order + 1) // 6, dtype=np_type)
    w = np.array([], dtype=np_type)

    coords = np.array([1.0, 0.0, 0.0,
                       2.0, 0.0, 0.0,
                       0.0, 1.0, 0.0,
                       0.0, 0.0, 1.0], dtype=np.float64)
    default_integral.tabulate_tensor(
        ffi.cast('{type} *'.format(type=c_type), b.ctypes.data),
        ffi.cast('{type} *'.format(type=c_type), w.ctypes.data),
        ffi.NULL,
        ffi.cast('double *', coords.ctypes.data), ffi.NULL, ffi.NULL, 0)

    # Check that the result is the same as for sympy
    assert np.allclose(b, [float(i) for i in sym])

    # Check that passing in permutations correctly reverses edges
    for edge in range(6):
        perm_b = np.zeros((order + 3) * (order + 2) * (order + 1) // 6, dtype=np_type)
        default_integral.tabulate_tensor(
            ffi.cast('{type} *'.format(type=c_type), perm_b.ctypes.data),
            ffi.cast('{type} *'.format(type=c_type), w.ctypes.data),
            ffi.NULL,
            ffi.cast('double *', coords.ctypes.data), ffi.NULL, ffi.NULL, 1 << (12 + edge))

        for e in range(6):
            start = 4 + (order - 1) * e
            end = start + order - 1
            if e == edge:
                assert np.allclose(b[start: end], perm_b[end - 1: start - 1: -1])
            else:
                assert np.allclose(b[start: end], perm_b[start: end])

    # Check that passing in permutation conts correctly permutes face data
    start = 4 + 6 * (order - 1)
    end = start + (order - 1) * (order - 2) // 2
    for rots in range(3):
        for refs in range(2):
            perm = (rots * 2 + refs)  # << 6
            new_coords = [1., 0., 0.]
            points = [[2.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
            new_coords += points[rots]
            if refs:
                new_coords += points[(rots - 1) % 3]
                new_coords += points[(rots + 1) % 3]
            else:
                new_coords += points[(rots + 1) % 3]
                new_coords += points[(rots - 1) % 3]
            new_coords = np.array(new_coords, dtype=np.float64)
            perm_b = np.zeros((order + 3) * (order + 2) * (order + 1) // 6, dtype=np_type)
            default_integral.tabulate_tensor(
                ffi.cast('{type} *'.format(type=c_type), perm_b.ctypes.data),
                ffi.cast('{type} *'.format(type=c_type), w.ctypes.data),
                ffi.NULL,
                ffi.cast('double *', new_coords.ctypes.data), ffi.NULL, ffi.NULL, perm)
            assert np.allclose(b[start:end], perm_b[start:end])
