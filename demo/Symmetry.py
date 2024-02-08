"""Symmetry demo."""

import basix.ufl
from ufl import FunctionSpace, Mesh, TestFunction, TrialFunction, dx, grad, inner

P1 = basix.ufl.element("P", "triangle", 1, shape=(2, 2), symmetry=True)
domain = Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
space = FunctionSpace(domain, P1)

u = TrialFunction(space)
v = TestFunction(space)

a = inner(grad(u), grad(v)) * dx
