# Copyright (C) 2016 Lizao Li
"""Biharmonis HHJ demo.

The bilinear form a(u, v) and linear form L(v) for Biharmonic equation
in Hellan-Herrmann-Johnson (HHJ) formulation.
"""

import basix.ufl
from ufl import (
    Coefficient,
    FacetNormal,
    FunctionSpace,
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
)

HHJ = basix.ufl.element("HHJ", "triangle", 2)
P = basix.ufl.element("P", "triangle", 3)
mixed_element = basix.ufl.mixed_element([HHJ, P])
domain = Mesh(basix.ufl.element("P", "triangle", 1, shape=(2,)))
mixed_space = FunctionSpace(domain, mixed_element)
p_space = FunctionSpace(domain, P)

(sigma, u) = TrialFunctions(mixed_space)
(tau, v) = TestFunctions(mixed_space)
f = Coefficient(p_space)


def b(sigma, v):
    """The form b."""
    n = FacetNormal(domain)
    return (
        inner(sigma, grad(grad(v))) * dx
        - dot(dot(sigma("+"), n("+")), n("+")) * jump(grad(v), n) * dS
        - dot(dot(sigma, n), n) * dot(grad(v), n) * ds
    )


a = inner(sigma, tau) * dx - b(tau, u) + b(sigma, v)
L = f * v * dx
