from ufl import (FiniteElement, TestFunction, TrialFunction, dx, inner,
                 quadrilateral)

element = FiniteElement("Lagrange", quadrilateral, 1)

v = TestFunction(element)
u = TrialFunction(element)

a = inner(u, v) * dx
