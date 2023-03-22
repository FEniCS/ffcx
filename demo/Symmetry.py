from basix.ufl import element
from ufl import TestFunction, TrialFunction, dx, grad, inner

P1 = element("P", "triangle", 1, shape=(2, 2), symmetry=True)

u = TrialFunction(P1)
v = TestFunction(P1)

a = inner(grad(u), grad(v)) * dx
