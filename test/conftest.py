# Copyright (C) 2020 Michal Habera
#
# This file is part of FFCX.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import pytest

import ffcx
import ufl

elements = [("Lagrange", ufl.interval, 1),
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
            # ("Lagrange", ufl.quadrilateral, 3, None, None, "spectral"),
            ("Lagrange", ufl.hexahedron, 3),
            # ("Lagrange", ufl.hexahedron, 3, None, None, "spectral"),
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
            ("N2curl", ufl.tetrahedron, 3),
            # ("Quadrature", ufl.interval, 2, None, "default"),
            # ("Quadrature", ufl.triangle, 2, None, "default"),
            # ("Quadrature", ufl.tetrahedron, 2, None, "default"),
            # ("Quadrature", ufl.quadrilateral, 2, None, "default"),
            # ("Quadrature", ufl.hexahedron, 2, None, "default")
            ]


@pytest.fixture(scope="module")
def compile_args():
    return ["-O0", "-Wall", "-Werror"]


@pytest.fixture(scope="module", params=elements, ids=["{!s}".format(el) for el in elements])
def compiled_element(compile_args, request):
    """Precompiled finite elements.

    Returns
    -------
    {(family, cell, degree): (ufl_element, compiled_element, compiled_module)}

    """
    element = request.param

    ufl_element = ufl.FiniteElement(*element)
    jit_compiled_elements, module = ffcx.codegeneration.jit.compile_elements(
        [ufl_element], cffi_extra_compile_args=compile_args)
    return (ufl_element, jit_compiled_elements[0], module)
