#
# Author: Martin Sandve Alnes
# Date: 2008-12-22
#

import basix.ufl
# Modified by Garth N. Wells, 2009
from ufl import (Coefficient, Constant, FacetNormal, Identity,
                 SpatialCoordinate, TestFunction, TrialFunction, derivative,
                 det, diff, dot, ds, dx, exp, grad, inner, inv, tetrahedron,
                 tr, variable)

# Cell and its properties
cell = tetrahedron
d = cell.geometric_dimension()
N = FacetNormal(cell)
x = SpatialCoordinate(cell)

# Elements
u_element = basix.ufl.element("P", cell.cellname(), 2, shape=(3, ))
p_element = basix.ufl.element("P", cell.cellname(), 1)
A_element = basix.ufl.element("P", cell.cellname(), 1, shape=(3, 3))

# Test and trial functions
v = TestFunction(u_element)
w = TrialFunction(u_element)

# Displacement at current and two previous timesteps
u = Coefficient(u_element)
up = Coefficient(u_element)
upp = Coefficient(u_element)

# Time parameters
dt = Constant(cell)

# Fiber field
A = Coefficient(A_element)

# External forces
T = Coefficient(u_element)
p0 = Coefficient(p_element)

# Material parameters FIXME
rho = Constant(cell)
K = Constant(cell)
c00 = Constant(cell)
c11 = Constant(cell)
c22 = Constant(cell)

# Deformation gradient
I = Identity(d)
F = I + grad(u)
F = variable(F)
Finv = inv(F)
J = det(F)

# Left Cauchy-Green deformation tensor
B = F * F.T
I1_B = tr(B)
I2_B = (I1_B**2 - tr(B * B)) / 2
I3_B = J**2

# Right Cauchy-Green deformation tensor
C = F.T * F
I1_C = tr(C)
I2_C = (I1_C**2 - tr(C * C)) / 2
I3_C = J**2

# Green strain tensor
E = (C - I) / 2

# Mapping of strain in fiber directions
Ef = A * E * A.T

# Strain energy function W(Q(Ef))
Q = c00 * Ef[0, 0]**2 + c11 * Ef[1, 1]**2 + c22 * Ef[2, 2]**2  # FIXME: insert some simple law here
W = (K / 2) * (exp(Q) - 1)  # + p stuff

# First Piola-Kirchoff stress tensor
P = diff(W, F)

# Acceleration term discretized with finite differences
k = dt / rho
acc = (u - 2 * up + upp)

# Residual equation # FIXME: Can contain errors, not tested!
a_F = inner(acc, v) * dx \
    + k * inner(P, grad(v)) * dx \
    - k * dot(J * Finv * T, v) * ds(0) \
    - k * dot(J * Finv * p0 * N, v) * ds(1)

# Jacobi matrix of residual equation
a_J = derivative(a_F, u, w)

# Export forms
forms = [a_F, a_J]
