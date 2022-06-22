from ufl import (FiniteElement, TestFunction, Coefficient, dx, inner,
                 quadrilateral)

element = FiniteElement("Lagrange", quadrilateral, 1)

v = TestFunction(element)
u = Coefficient(element)

a = inner(u, v) * dx
