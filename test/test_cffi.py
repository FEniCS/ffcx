# -*- coding: utf-8 -*-
# Copyright (C) 2018 Garth N. Wells
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import numpy as np
import pytest

import ffc.backends.ufc.jit
import ufl


@pytest.fixture(scope="module")
def lagrange_element():
    """Compile list of Lagrange elements"""
    cell = ufl.triangle
    elements = [ufl.FiniteElement("Lagrange", cell, p) for p in range(1, 5)]
    compiled_elements, module = ffc.backends.ufc.jit.compile_elements(elements)
    return elements, compiled_elements, module


def test_dim_degree(lagrange_element):
    ufl_elements, compiled_elements, module = lagrange_element
    for e, compiled_e in zip(ufl_elements, compiled_elements):
        assert compiled_e.geometric_dimension == 2
        assert compiled_e.topological_dimension == 2
        assert e.degree() == compiled_e.degree


def test_tabulate_reference_dof_coordinates(lagrange_element):
    ufl_elements, compiled_elements, module = lagrange_element
    for e, compiled_e in zip(ufl_elements, compiled_elements):
        #test = module.ffi.string(compiled_e.family)
        tdim = compiled_e.topological_dimension
        space_dim = compiled_e.space_dimension
        X = np.zeros([space_dim, tdim])
        X_ptr = module.ffi.cast("double *", module.ffi.from_buffer(X))
        compiled_e.tabulate_reference_dof_coordinates(X_ptr)
        # print(X)


def test_form():
    cell = ufl.triangle
    element = ufl.FiniteElement("Lagrange", cell, 1)
    u, v = ufl.TestFunction(element), ufl.TrialFunction(element),
    a = ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx
    forms = [a]
    compiled_forms, module = ffc.backends.ufc.jit.compile_forms(forms)

    for f, compiled_f in zip(forms, compiled_forms):
        assert compiled_f.rank == f.rank()


# cell = ufl.triangle
# elements = [ufl.FiniteElement("Lagrange", cell, p) for p in range(1, 5)]
# compiled_elements, module = ffc.backends.ufc.jit.compile_elements(elements)

# for e, compiled_e in zip(elements, compiled_elements):
#     assert compiled_e.geometric_dimension == 2
#     assert compiled_e.topological_dimension == 2
#     assert e.degree() == compiled_e.degree
