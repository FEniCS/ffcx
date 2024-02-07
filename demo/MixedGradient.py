"""Mixed gradient demo."""

import basix.ufl
from ufl import FunctionSpace, Mesh, TestFunctions, TrialFunctions, ds, grad, inner

element1 = basix.ufl.element("DG", "triangle", 1)
element2 = basix.ufl.element("DGT", "triangle", 1)
element = basix.ufl.mixed_element([element1, element2])
domain = Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
space = FunctionSpace(domain, element)

u = TrialFunctions(space)[0]
v = TestFunctions(space)[0]

a = inner(grad(u), grad(v)) * ds
