from ufl import (FiniteElement, TestFunction, Coefficient, TrialFunction, dx, inner,
                 hexahedron)

element = FiniteElement("Lagrange", hexahedron, 5)

v = TestFunction(element)
u = TrialFunction(element)

a = inner(u, v) * dx
