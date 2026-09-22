# Copyright (C) 2026 Garth N. Wells
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Stable generated kernel name demo.

The bilinear form uses the Python variable name to derive its kernel name.
The linear form demonstrates an explicit ``ffcx_kernel_name`` override.
"""

import basix.ufl
from ufl import FunctionSpace, Mesh, TestFunction, TrialFunction, dx, inner

element = basix.ufl.element("Lagrange", "triangle", 1)
domain = Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
space = FunctionSpace(domain, element)

u = TrialFunction(space)
v = TestFunction(space)

mass = inner(u, v) * dx
unit_load = inner(1.0, v) * dx(metadata={"ffcx_kernel_name": "tabulate_tensor_unit_load"})

forms = [mass, unit_load]
