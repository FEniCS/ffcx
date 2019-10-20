# Copyright (C) 2019 Chris Richardson
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import numpy as np
import pytest
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


@pytest.mark.parametrize("mode", ["double", "float", "long double", "double complex", "float complex"])
def test_additive_facet_integral(mode):
    cell = ufl.triangle
    element = ufl.FiniteElement("Lagrange", cell, 1)
    u, v = ufl.TrialFunction(element), ufl.TestFunction(element)
    a = ufl.inner(u, v) * ufl.ds
    forms = [a]
    compiled_forms, module = ffc.codegeneration.jit.compile_forms(forms, parameters={'scalar_type': mode})

    for f, compiled_f in zip(forms, compiled_forms):
        assert compiled_f.rank == len(f.arguments())

    ffi = cffi.FFI()
    form0 = compiled_forms[0][0]

    assert form0.num_exterior_facet_integrals == 1
    ids = np.zeros(form0.num_exterior_facet_integrals, dtype=np.int32)
    form0.get_exterior_facet_integral_ids(ffi.cast('int *', ids.ctypes.data))
    assert ids[0] == -1

    default_integral = form0.create_exterior_facet_integral(ids[0])

    c_type, np_type = float_to_type(mode)
    A = np.zeros((3, 3), dtype=np_type)
    w = np.array([], dtype=np_type)
    c = np.array([], dtype=np_type)
    facets = np.array([0], dtype=np.int32)

    coords = np.array([0.0, 2.0, np.sqrt(3.0), -1.0, -np.sqrt(3.0), -1.0], dtype=np.float64)

    for i in range(3):
        facets[0] = i
        default_integral.tabulate_tensor(
            ffi.cast('{type} *'.format(type=c_type), A.ctypes.data),
            ffi.cast('{type} *'.format(type=c_type), w.ctypes.data),
            ffi.cast('{type} *'.format(type=c_type), c.ctypes.data),
            ffi.cast('double *', coords.ctypes.data),
            ffi.cast('int *', facets.ctypes.data), ffi.NULL)

        assert np.isclose(A.sum(), np.sqrt(12) * (i + 1))


@pytest.mark.parametrize("mode", ["double", "float", "long double", "double complex", "float complex"])
def test_additive_cell_integral(mode):
    cell = ufl.triangle
    element = ufl.FiniteElement("Lagrange", cell, 1)
    u, v = ufl.TrialFunction(element), ufl.TestFunction(element)
    a = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx
    forms = [a]
    compiled_forms, module = ffc.codegeneration.jit.compile_forms(forms, parameters={'scalar_type': mode})

    for f, compiled_f in zip(forms, compiled_forms):
        assert compiled_f.rank == len(f.arguments())

    ffi = cffi.FFI()
    form0 = compiled_forms[0][0]

    assert form0.num_cell_integrals == 1
    ids = np.zeros(form0.num_cell_integrals, dtype=np.int32)
    form0.get_cell_integral_ids(ffi.cast('int *', ids.ctypes.data))
    assert ids[0] == -1

    default_integral = form0.create_cell_integral(ids[0])

    c_type, np_type = float_to_type(mode)
    A = np.zeros((3, 3), dtype=np_type)
    w = np.array([], dtype=np_type)
    c = np.array([], dtype=np_type)

    coords = np.array([0.0, 2.0, np.sqrt(3.0), -1.0, -np.sqrt(3.0), -1.0], dtype=np.float64)

    default_integral.tabulate_tensor(
        ffi.cast('{type} *'.format(type=c_type), A.ctypes.data),
        ffi.cast('{type} *'.format(type=c_type), w.ctypes.data),
        ffi.cast('{type} *'.format(type=c_type), c.ctypes.data),
        ffi.cast('double *', coords.ctypes.data), ffi.NULL, ffi.NULL)

    A0 = np.array(A)

    for i in range(3):

        default_integral.tabulate_tensor(
            ffi.cast('{type} *'.format(type=c_type), A.ctypes.data),
            ffi.cast('{type} *'.format(type=c_type), w.ctypes.data),
            ffi.cast('{type} *'.format(type=c_type), c.ctypes.data),
            ffi.cast('double *', coords.ctypes.data), ffi.NULL, ffi.NULL)

        assert np.all(np.isclose(A, (i + 2) * A0))
