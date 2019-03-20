# -*- coding: utf-8 -*-
# Copyright (C) 2019 Chris Richardson
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import numpy as np
import cffi

import ffc.codegeneration.jit
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


def test_custom_form_2d():
    mode = 'double'
    cell = ufl.triangle
    element = ufl.FiniteElement("Lagrange", cell, 1)
    compiled_element, module = ffc.codegeneration.jit.compile_elements([element])
    el = compiled_element[0][0]
    sd = el.space_dimension

    u, v = ufl.TrialFunction(element), ufl.TestFunction(element)
    a1 = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dc + ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx
    a2 = u * v * ufl.dc + u * v * ufl.dx
    forms = [a1, a2]
    compiled_forms, module = ffc.codegeneration.jit.compile_forms(
        forms, parameters={'scalar_type': mode})

    for f, compiled_f in zip(forms, compiled_forms):
        assert compiled_f.rank == len(f.arguments())

    for f in compiled_forms:
        form0 = f[0].create_custom_integral(-1)
        form1 = f[0].create_cell_integral(-1)

        c_type, np_type = float_to_type(mode)
        A0 = np.zeros((sd, sd), dtype=np_type)
        w = np.array([], dtype=np_type)
        ffi = cffi.FFI()
        coords = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]], dtype=np.float64)

        # Quadrature points and weights
        Qpts = np.array([[0.81684757, 0.09157621],
                         [0.09157621, 0.81684757],
                         [0.09157621, 0.09157621],
                         [0.10810302, 0.44594849],
                         [0.44594849, 0.10810302],
                         [0.44594849, 0.44594849]])
        Qwts = np.array([0.05497587, 0.05497587, 0.05497587, 0.11169079, 0.11169079, 0.11169079])

        Fnormal = np.array([])

        form0.tabulate_tensor(
            ffi.cast('{type} *'.format(type=c_type), A0.ctypes.data),
            ffi.cast('{type} *'.format(type=c_type), w.ctypes.data),
            ffi.cast('double *', coords.ctypes.data), len(Qpts),
            ffi.cast('{type} *'.format(type=c_type), Qpts.ctypes.data),
            ffi.cast('{type} *'.format(type=c_type), Qwts.ctypes.data),
            ffi.cast('{type} *'.format(type=c_type), Fnormal.ctypes.data), 0)

        A1 = np.zeros((sd, sd), dtype=np_type)

        form1.tabulate_tensor(
            ffi.cast('{type} *'.format(type=c_type), A1.ctypes.data),
            ffi.cast('{type} *'.format(type=c_type), w.ctypes.data),
            ffi.cast('double *', coords.ctypes.data), 0)

        print(A0, A1)

        assert np.isclose(A0, A1).all()
