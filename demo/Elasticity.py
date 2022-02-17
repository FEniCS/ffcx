#
# Author: Anders Logg
# Modified by: Martin Sandve Alnes
# Date: 2009-01-12
#
from ufl import (TestFunction, TrialFunction, VectorElement, dx, grad, inner,
                 tetrahedron)

element = VectorElement("Lagrange", tetrahedron, 1)

v = TestFunction(element)
u = TrialFunction(element)


def epsilon(v):
    Dv = grad(v)
    return 0.5 * (Dv + Dv.T)


a = inner(epsilon(v), epsilon(u)) * dx
