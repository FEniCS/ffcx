# Copyright (C) 2016 Lizao Li
"""Biharmonic Regge demo.

The bilinear form a(u, v) and linear form L(v) for Biharmonic equation in Regge formulation.
"""

import basix.ufl
from ufl import (
    Coefficient,
    FacetNormal,
    FunctionSpace,
    Identity,
    Mesh,
    TestFunctions,
    TrialFunctions,
    dot,
    dS,
    ds,
    dx,
    grad,
    inner,
    jump,
    tr,
)

REG = basix.ufl.element("Regge", "tetrahedron", 1)
P = basix.ufl.element("Lagrange", "tetrahedron", 2)
mixed_element = basix.ufl.mixed_element([REG, P])
domain = Mesh(basix.ufl.element("P", "tetrahedron", 1, shape=(3,)))
mixed_space = FunctionSpace(domain, mixed_element)
p_space = FunctionSpace(domain, P)

(sigma, u) = TrialFunctions(mixed_space)
(tau, v) = TestFunctions(mixed_space)
f = Coefficient(p_space)


def S(mu):
    """The form S."""
    return mu - Identity(3) * tr(mu)


def b(mu, v):
    """The form b."""
    n = FacetNormal(domain)
    return (
        inner(S(mu), grad(grad(v))) * dx
        - dot(dot(S(mu("+")), n("+")), n("+")) * jump(grad(v), n) * dS
        - dot(dot(S(mu), n), n) * dot(grad(v), n) * ds
    )


a = inner(S(sigma), S(tau)) * dx - b(tau, u) + b(sigma, v)
L = f * v * dx
