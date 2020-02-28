# Copyright (C) 2020 Michal Habera
#
# This file is part of FFCX.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import pytest

import ufl
import ffcx


@pytest.fixture(scope="module")
def elements():
    return [("Lagrange", ufl.interval, 1),
            ("Lagrange", ufl.triangle, 1),
            ("Lagrange", ufl.tetrahedron, 1),
            ("Lagrange", ufl.quadrilateral, 1),
            ("Lagrange", ufl.hexahedron, 1),
            ("Lagrange", ufl.interval, 2),
            ("Lagrange", ufl.triangle, 2),
            ("Lagrange", ufl.tetrahedron, 2),
            ("Lagrange", ufl.quadrilateral, 2),
            ("Lagrange", ufl.hexahedron, 2),
            ("Lagrange", ufl.interval, 3),
            ("Lagrange", ufl.triangle, 3),
            ("Lagrange", ufl.tetrahedron, 3),
            ("Lagrange", ufl.quadrilateral, 3),
            ("Lagrange", ufl.hexahedron, 3),
            ("Brezzi-Douglas-Marini", ufl.triangle, 1),
            ("Brezzi-Douglas-Marini", ufl.tetrahedron, 1),
            ("Brezzi-Douglas-Marini", ufl.triangle, 2),
            ("Brezzi-Douglas-Marini", ufl.tetrahedron, 2),
            ("Brezzi-Douglas-Marini", ufl.triangle, 3),
            ("Brezzi-Douglas-Marini", ufl.tetrahedron, 3),
            ("Raviart-Thomas", ufl.triangle, 1),
            ("Raviart-Thomas", ufl.tetrahedron, 1),
            ("Raviart-Thomas", ufl.triangle, 2),
            ("Raviart-Thomas", ufl.tetrahedron, 2),
            ("Raviart-Thomas", ufl.triangle, 3),
            ("Raviart-Thomas", ufl.tetrahedron, 3),
            ("N1curl", ufl.triangle, 1),
            ("N1curl", ufl.tetrahedron, 1),
            ("N1curl", ufl.triangle, 2),
            ("N1curl", ufl.tetrahedron, 2),
            ("N1curl", ufl.triangle, 3),
            ("N1curl", ufl.tetrahedron, 3),
            ("N2curl", ufl.triangle, 1),
            ("N2curl", ufl.tetrahedron, 1),
            ("N2curl", ufl.triangle, 2),
            ("N2curl", ufl.tetrahedron, 2),
            ("N2curl", ufl.triangle, 3),
            ("N2curl", ufl.tetrahedron, 3)]


@pytest.fixture(scope="module")
def compile_args():
    return ["-O0", "-Wall", "-Werror"]


@pytest.fixture(scope="module")
def compiled_elements(compile_args, elements):
    comp_elements = {}
    for element in elements:
        ufl_element = ufl.FiniteElement(*element)
        jit_compiled_elements, module = ffcx.codegeneration.jit.compile_elements(
            [ufl_element], cffi_extra_compile_args=compile_args)
        comp_elements[element] = (ufl_element, jit_compiled_elements[0], module)

    return comp_elements
