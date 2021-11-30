# Copyright (C) 2021 Igor A. Baratta
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later


from ffcx.codegeneration.flop_count import count_flops
import ufl


def create_form(degree):
    mesh = ufl.Mesh(ufl.VectorElement("Lagrange", "triangle", 1))
    element = ufl.FiniteElement("Lagrange", ufl.triangle, degree)
    V = ufl.FunctionSpace(mesh, element)

    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)

    return ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx + ufl.inner(u, v) * ufl.ds


def test_flops():
    k1, k2 = 2, 4
    a1 = create_form(k1)
    a2 = create_form(k2)

    dofs1 = (k1 + 1.) * (k1 + 2.) / 2.
    dofs2 = (k2 + 1.) * (k2 + 2.) / 2.

    flops_1 = count_flops(a1)
    assert(len(flops_1) == 2)

    flops_2 = count_flops(a2)
    assert(len(flops_2) == 2)

    r = sum(flops_2, 0.) / sum(flops_1, 0.)

    assert r > (dofs2**2 / dofs1**2)
