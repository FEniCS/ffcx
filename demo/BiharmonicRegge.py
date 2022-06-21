# Copyright (C) 2016 Lizao Li
#
# The bilinear form a(u, v) and linear form L(v) for
# Biharmonic equation in Regge formulation.
from ufl import (Coefficient, FacetNormal, FiniteElement, Identity,
                 TestFunctions, TrialFunctions, dot, dS, ds, dx, grad, inner,
                 jump, tetrahedron, tr)

REG = FiniteElement('Regge', tetrahedron, 1)
CG = FiniteElement('Lagrange', tetrahedron, 2)
mixed_element = REG * CG

(sigma, u) = TrialFunctions(mixed_element)
(tau, v) = TestFunctions(mixed_element)
f = Coefficient(CG)


def S(mu):
    return mu - Identity(3) * tr(mu)


def b(mu, v):
    n = FacetNormal(tetrahedron)
    return inner(S(mu), grad(grad(v))) * dx \
        - dot(dot(S(mu('+')), n('+')), n('+')) * jump(grad(v), n) * dS \
        - dot(dot(S(mu), n), n) * dot(grad(v), n) * ds


a = inner(S(sigma), S(tau)) * dx - b(tau, u) + b(sigma, v)
L = f * v * dx
