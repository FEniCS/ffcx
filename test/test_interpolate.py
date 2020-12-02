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


@pytest.mark.parametrize("mode", ["double"])
@pytest.mark.parametrize("order", [1, 2, 3, 4, 5])
def test_lagrange_triangle(compile_args, order, mode):
    cell = ufl.triangle
    element = ufl.FiniteElement("Lagrange", cell, order)
    compiled_elements, module = ffcx.codegeneration.jit.compile_elements(
        [element], parameters={'scalar_type': mode}, cffi_extra_compile_args=compile_args)

    ffi = cffi.FFI()
    element0 = compiled_elements[0][0]

    c_type, np_type = float_to_type(mode)
    coeffs_in = np.array(range(element0.space_dimension), dtype=np_type)

    coeffs = np.zeros((order + 2) * (order + 1) // 2, dtype=np_type)

    coords = np.array([1.0, 0.0, 2.0, 0.0, 0.0, 1.0], dtype=np.float64)
    element0.interpolate_into_cell(
        ffi.cast('{type} *'.format(type=c_type), coeffs.ctypes.data),
        ffi.cast('{type} *'.format(type=c_type), coeffs_in.ctypes.data), 0)

    # Check that the result is the same as the input
    assert np.allclose(coeffs, coeffs_in)

    # Check that passing in permutations correctly reverses edges
    for perm in range(1, 8):
        element0.interpolate_into_cell(
            ffi.cast('{type} *'.format(type=c_type), coeffs.ctypes.data),
            ffi.cast('{type} *'.format(type=c_type), coeffs_in.ctypes.data),
            perm)

        for edge in range(3):
            start = 3 + (order - 1) * edge
            end = start + order - 1
            if perm >> edge & 1 == 1:
                assert np.allclose(coeffs[start: end], range(end - 1, start - 1, -1))
            else:
                assert np.allclose(coeffs[start: end], range(start, end))


@pytest.mark.parametrize("mode", ["double"])
@pytest.mark.parametrize("order", [1, 2, 3, 4])
def test_lagrange_tetrahedron(compile_args, order, mode):
    cell = ufl.tetrahedron
    element = ufl.FiniteElement("Lagrange", cell, order)
    compiled_elements, module = ffcx.codegeneration.jit.compile_elements(
        [element], parameters={'scalar_type': mode}, cffi_extra_compile_args=compile_args)

    ffi = cffi.FFI()
    element0 = compiled_elements[0][0]

    c_type, np_type = float_to_type(mode)
    coeffs_in = np.array(range(element0.space_dimension), dtype=np_type)

    coeffs = np.zeros((order + 3) * (order + 2) * (order + 1) // 6, dtype=np_type)

    coords = np.array([1.0, 0.0, 2.0, 0.0, 0.0, 1.0], dtype=np.float64)
    element0.interpolate_into_cell(
        ffi.cast('{type} *'.format(type=c_type), coeffs.ctypes.data),
        ffi.cast('{type} *'.format(type=c_type), coeffs_in.ctypes.data), 0)

    # Check that the result is the same as the input
    assert np.allclose(coeffs, coeffs_in)

    # TODO: debug on TETS
    # Check that passing in permutations correctly reverses edges
    for edge in range(6):
        element0.interpolate_into_cell(
            ffi.cast('{type} *'.format(type=c_type), coeffs.ctypes.data),
            ffi.cast('{type} *'.format(type=c_type), coeffs_in.ctypes.data),
            1 << (12 + edge))

        for e in range(6):
            start = 4 + (order - 1) * e
            end = start + order - 1
            if e == edge:
                assert np.allclose(coeffs[start: end], range(end - 1, start - 1, -1))
            else:
                assert np.allclose(coeffs[start: end], range(start, end))


    # Check that passing in permutation conts correctly permutes face data
    if order > 4:
        raise NotImplementedError("This test is not implemented yet for higher of order >4.")
    if order == 4:
        expected = [(0, 1, 2), (0, 2, 1),
                    (2, 0, 1), (2, 1, 0),
                    (1, 2, 0), (1, 0, 2)]
        for face in range(4):
            start = 4 + 6 * (order - 1) + face * (order - 1) * (order - 2) // 2
            end = start + (order - 1) * (order - 2) // 2
            for rots in range(3):
                for refs in range(2):
                    perm = rots * 2 + refs
                    element0.interpolate_into_cell(
                        ffi.cast('{type} *'.format(type=c_type), coeffs.ctypes.data),
                        ffi.cast('{type} *'.format(type=c_type), coeffs_in.ctypes.data),
                        perm << (3 * face))

                    assert np.allclose(coeffs[:start], range(start))
                    assert np.allclose(coeffs[start:end], [start + i for i in expected[perm]])
                    assert np.allclose(coeffs[end:], range(end, len(coeffs)))
