from ufl import (TensorElement, TestFunction, TrialFunction, dx, grad, inner,
                 triangle)

P1 = TensorElement("CG", triangle, 1, (2, 2), symmetry={(1, 0): (0, 1)})

u = TrialFunction(P1)
v = TestFunction(P1)

a = inner(grad(u), grad(v)) * dx
