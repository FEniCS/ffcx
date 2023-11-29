# Copyright (C) 2019 Chris Richardson
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import numpy as np
import pytest

import ffcx.codegeneration.jit
import basix.ufl
import ufl
from ffcx.codegeneration.utils import cdtype_to_numpy, scalar_to_value_type


@pytest.mark.parametrize("mode",
                         [
                             "double",
                             "float",
                             "long double",
                             "double _Complex",
                             "float _Complex"
                         ])
def test_additive_facet_integral(mode, compile_args):
    element = basix.ufl.element("Lagrange", "triangle", 1)
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2, )))
    space = ufl.FunctionSpace(domain, element)
    u, v = ufl.TrialFunction(space), ufl.TestFunction(space)
    a = ufl.inner(u, v) * ufl.ds
    forms = [a]
    compiled_forms, module, code = ffcx.codegeneration.jit.compile_forms(
        forms, options={'scalar_type': mode}, cffi_extra_compile_args=compile_args)

    for f, compiled_f in zip(forms, compiled_forms):
        assert compiled_f.rank == len(f.arguments())

    ffi = module.ffi
    form0 = compiled_forms[0]

    integral_offsets = form0.form_integral_offsets
    ex = module.lib.exterior_facet
    assert integral_offsets[ex + 1] - integral_offsets[ex] == 1
    integral_id = form0.form_integral_ids[integral_offsets[ex]]
    assert integral_id == -1

    default_integral = form0.form_integrals[integral_offsets[ex]]

    np_type = cdtype_to_numpy(mode)
    A = np.zeros((3, 3), dtype=np_type)
    w = np.array([], dtype=np_type)
    c = np.array([], dtype=np_type)
    facets = np.array([0], dtype=np.int32)
    perm = np.array([0], dtype=np.uint8)

    geom_type = scalar_to_value_type(mode)
    np_gtype = cdtype_to_numpy(geom_type)
    coords = np.array([0.0, 2.0, 0.0,
                       np.sqrt(3.0), -1.0, 0.0,
                       -np.sqrt(3.0), -1.0, 0.0], dtype=np_gtype)

    kernel = getattr(default_integral, f"tabulate_tensor_{np_type}")

    for i in range(3):
        facets[0] = i
        kernel(ffi.cast('{type} *'.format(type=mode), A.ctypes.data),
               ffi.cast('{type} *'.format(type=mode), w.ctypes.data),
               ffi.cast('{type} *'.format(type=mode), c.ctypes.data),
               ffi.cast(f'{geom_type} *', coords.ctypes.data),
               ffi.cast('int *', facets.ctypes.data),
               ffi.cast('uint8_t *', perm.ctypes.data))

        assert np.isclose(A.sum(), np.sqrt(12) * (i + 1))


@pytest.mark.parametrize("mode", ["double", "float", "long double", "double _Complex", "float _Complex"])
def test_additive_cell_integral(mode, compile_args):
    element = basix.ufl.element("Lagrange", "triangle", 1)
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2, )))
    space = ufl.FunctionSpace(domain, element)
    u, v = ufl.TrialFunction(space), ufl.TestFunction(space)
    a = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx
    forms = [a]
    compiled_forms, module, code = ffcx.codegeneration.jit.compile_forms(
        forms, options={'scalar_type': mode}, cffi_extra_compile_args=compile_args)

    for f, compiled_f in zip(forms, compiled_forms):
        assert compiled_f.rank == len(f.arguments())

    ffi = module.ffi
    form0 = compiled_forms[0]

    cell = module.lib.cell
    offsets = form0.form_integral_offsets
    num_integrals = offsets[cell + 1] - offsets[cell]
    assert num_integrals == 1
    integral_id = form0.form_integral_ids[offsets[cell]]
    assert integral_id == -1

    default_integral = form0.form_integrals[offsets[cell]]

    np_type = cdtype_to_numpy(mode)
    A = np.zeros((3, 3), dtype=np_type)
    w = np.array([], dtype=np_type)
    c = np.array([], dtype=np_type)

    geom_type = scalar_to_value_type(mode)
    np_gtype = cdtype_to_numpy(geom_type)
    coords = np.array([0.0, 2.0, 0.0,
                       np.sqrt(3.0), -1.0, 0.0,
                       -np.sqrt(3.0), -1.0, 0.0], dtype=np_gtype)

    kernel = getattr(default_integral, f"tabulate_tensor_{np_type}")

    kernel(ffi.cast('{type} *'.format(type=mode), A.ctypes.data),
           ffi.cast('{type} *'.format(type=mode), w.ctypes.data),
           ffi.cast('{type} *'.format(type=mode), c.ctypes.data),
           ffi.cast(f'{geom_type} *', coords.ctypes.data), ffi.NULL, ffi.NULL)

    A0 = np.array(A)
    for i in range(3):
        kernel(ffi.cast('{type} *'.format(type=mode), A.ctypes.data),
               ffi.cast('{type} *'.format(type=mode), w.ctypes.data),
               ffi.cast('{type} *'.format(type=mode), c.ctypes.data),
               ffi.cast(f'{geom_type} *', coords.ctypes.data), ffi.NULL, ffi.NULL)

        assert np.all(np.isclose(A, (i + 2) * A0))
