# Copyright (C) 2026 Garth N. Wells
"""Poisson demo for the C example.

The bilinear form a(u, v) and linear form L(v) for Poisson's equation,
compiled by the CMake project in this directory.
"""

import basix.ufl
from ufl import Coefficient, FunctionSpace, Mesh, TestFunction, TrialFunction, dx, grad, inner

e = basix.ufl.element("Lagrange", "triangle", 1)
coord_element = basix.ufl.element("Lagrange", "triangle", 1, shape=(2,))
mesh = Mesh(coord_element)
V = FunctionSpace(mesh, e)

u = TrialFunction(V)
v = TestFunction(V)
f = Coefficient(V)

a = inner(grad(u), grad(v)) * dx
L = inner(f, v) * dx
forms = [a, L]
