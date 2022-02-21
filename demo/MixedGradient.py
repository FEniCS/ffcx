from ufl import (FiniteElement, MixedElement, TestFunctions, TrialFunctions,
                 ds, grad, inner, triangle)

element1 = FiniteElement("DG", triangle, 1)
element2 = FiniteElement("DGT", triangle, 1)
element = MixedElement(element1, element2)

u = TrialFunctions(element)[0]
v = TestFunctions(element)[0]

a = inner(grad(u), grad(v)) * ds
