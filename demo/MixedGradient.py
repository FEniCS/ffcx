import basix.ufl
from ufl import TestFunctions, TrialFunctions, ds, grad, inner

element1 = basix.ufl.element("DG", "triangle", 1)
element2 = basix.ufl.element("DGT", "triangle", 1)
element = basix.ufl.mixed_element([element1, element2])

u = TrialFunctions(element)[0]
v = TestFunctions(element)[0]

a = inner(grad(u), grad(v)) * ds
