
import basix.ufl
from ufl import (
    Coefficient,
    Constant,
    FunctionSpace,
    Mesh,
    TestFunction,
    TrialFunction,
    dx,
    grad,
    inner,
    tr,
)

mesh = Mesh(basix.ufl.element("P", "triangle", 1, shape=(2,)))

e = basix.ufl.element("Lagrange", "triangle", 1)
space = FunctionSpace(mesh, e)

u = TrialFunction(space)
v = TestFunction(space)
f = Coefficient(space)

kappa = Constant(mesh, shape=(2, 2))

a = tr(kappa) * inner(grad(u), grad(v)) * dx
# L = f * v * dx
