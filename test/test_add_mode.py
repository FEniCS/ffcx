# Copyright (C) 2019 Chris Richardson
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import numpy as np
import pytest

import ffcx.codegeneration.jit
import ufl
from ffcx.naming import cdtype_to_numpy


@pytest.mark.parametrize("mode", ["double", "float", "long double", "double _Complex", "float _Complex"])
def test_additive_facet_integral(mode, compile_args):
    cell = ufl.triangle
    element = ufl.FiniteElement("Lagrange", cell, 1)
    u, v = ufl.TrialFunction(element), ufl.TestFunction(element)
    a = ufl.inner(u, v) * ufl.ds
    forms = [a]
    compiled_forms, module, code = ffcx.codegeneration.jit.compile_forms(
        forms, parameters={'scalar_type': mode}, cffi_extra_compile_args=compile_args)

    for f, compiled_f in zip(forms, compiled_forms):
        assert compiled_f.rank == len(f.arguments())

    ffi = module.ffi
    form0 = compiled_forms[0]

    assert form0.num_integrals(module.lib.exterior_facet) == 1
    ids = form0.integral_ids(module.lib.exterior_facet)
    assert ids[0] == -1

    default_integral = form0.integrals(module.lib.exterior_facet)[0]

    np_type = cdtype_to_numpy(mode)
    A = np.zeros((3, 3), dtype=np_type)
    w = np.array([], dtype=np_type)
    c = np.array([], dtype=np_type)
    facets = np.array([0], dtype=np.int32)
    perm = np.array([0], dtype=np.uint8)

    coords = np.array([0.0, 2.0, 0.0,
                       np.sqrt(3.0), -1.0, 0.0,
                       -np.sqrt(3.0), -1.0, 0.0], dtype=np.float64)

    kernel = getattr(default_integral, f"tabulate_tensor_{np_type}")

    for i in range(3):
        facets[0] = i
        kernel(ffi.cast('{type} *'.format(type=mode), A.ctypes.data),
               ffi.cast('{type} *'.format(type=mode), w.ctypes.data),
               ffi.cast('{type} *'.format(type=mode), c.ctypes.data),
               ffi.cast('double *', coords.ctypes.data),
               ffi.cast('int *', facets.ctypes.data),
               ffi.cast('uint8_t *', perm.ctypes.data))

        assert np.isclose(A.sum(), np.sqrt(12) * (i + 1))


@pytest.mark.parametrize("mode", ["double", "float", "long double", "double _Complex", "float _Complex"])
def test_additive_cell_integral(mode, compile_args):
    cell = ufl.triangle
    element = ufl.FiniteElement("Lagrange", cell, 1)
    u, v = ufl.TrialFunction(element), ufl.TestFunction(element)
    a = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx
    forms = [a]
    compiled_forms, module, code = ffcx.codegeneration.jit.compile_forms(
        forms, parameters={'scalar_type': mode}, cffi_extra_compile_args=compile_args)

    for f, compiled_f in zip(forms, compiled_forms):
        assert compiled_f.rank == len(f.arguments())

    ffi = module.ffi
    form0 = compiled_forms[0]

    assert form0.num_integrals(module.lib.cell) == 1
    ids = form0.integral_ids(module.lib.cell)
    assert ids[0] == -1

    default_integral = form0.integrals(0)[0]

    np_type = cdtype_to_numpy(mode)
    A = np.zeros((3, 3), dtype=np_type)
    w = np.array([], dtype=np_type)
    c = np.array([], dtype=np_type)

    coords = np.array([0.0, 2.0, 0.0,
                       np.sqrt(3.0), -1.0, 0.0,
                       -np.sqrt(3.0), -1.0, 0.0], dtype=np.float64)

    kernel = getattr(default_integral, f"tabulate_tensor_{np_type}")

    kernel(ffi.cast('{type} *'.format(type=mode), A.ctypes.data),
           ffi.cast('{type} *'.format(type=mode), w.ctypes.data),
           ffi.cast('{type} *'.format(type=mode), c.ctypes.data),
           ffi.cast('double *', coords.ctypes.data), ffi.NULL, ffi.NULL)

    A0 = np.array(A)

    for i in range(3):
        kernel(ffi.cast('{type} *'.format(type=mode), A.ctypes.data),
               ffi.cast('{type} *'.format(type=mode), w.ctypes.data),
               ffi.cast('{type} *'.format(type=mode), c.ctypes.data),
               ffi.cast('double *', coords.ctypes.data), ffi.NULL, ffi.NULL)

        assert np.all(np.isclose(A, (i + 2) * A0))
