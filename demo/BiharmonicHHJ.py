# Copyright (C) 2016 Lizao Li
#
# The bilinear form a(u, v) and linear form L(v) for
# Biharmonic equation in Hellan-Herrmann-Johnson (HHJ)
# formulation.
import basix.ufl
from ufl import (Coefficient, FacetNormal, TestFunctions, TrialFunctions, dot,
                 dS, ds, dx, grad, inner, jump, triangle)

HHJ = basix.ufl.element('HHJ', "triangle", 2)
P = basix.ufl.element('P', "triangle", 3)
mixed_element = basix.ufl.mixed_element([HHJ, P])

(sigma, u) = TrialFunctions(mixed_element)
(tau, v) = TestFunctions(mixed_element)
f = Coefficient(P)


def b(sigma, v):
    n = FacetNormal(triangle)
    return inner(sigma, grad(grad(v))) * dx \
        - dot(dot(sigma('+'), n('+')), n('+')) * jump(grad(v), n) * dS \
        - dot(dot(sigma, n), n) * dot(grad(v), n) * ds


a = inner(sigma, tau) * dx - b(tau, u) + b(sigma, v)
L = f * v * dx
