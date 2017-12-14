# -*- coding: utf-8 -*-
"""
Tests of algorithm for factorization of integrand w.r.t. Argument terms.
"""

from ufl import *
from ffc.uflacs.analysis.factorization import compute_argument_factorization

# TODO: Restructure these tests using py.test fixtures and parameterization?


def compare_compute_argument_factorization(SV, SV_deps, expected_AV, expected_FV, expected_IM):
    SV_targets = [len(SV) - 1]
    rank = max(len(k) for k in expected_IM.keys())
    argument_factorizations, modified_arguments, FV, FV_deps, FV_targets = \
        compute_argument_factorization(SV, SV_deps, SV_targets, rank)
    argument_factorization, = argument_factorizations

    AV = modified_arguments
    IM = argument_factorization

    assert AV == expected_AV
    if 0:
        for n in range(1, min(len(FV), len(expected_FV))):
            assert FV[:n] == expected_FV[:n]
    assert FV == expected_FV
    assert IM == expected_IM


def test_compute_argument_factorization():
    V = FiniteElement("CG", triangle, 1)
    u = TrialFunction(V)
    v = TestFunction(V)
    a, b, c, d, e, f, g = [Coefficient(V, count=k) for k in range(7)]

    zero = as_ufl(0.0)
    one = as_ufl(1.0)
    two = as_ufl(2)

    # Test workaround for hack in factorization:
    FVpre = [zero, one, two]
    offset = len(FVpre)

    # Test basic non-argument terminal
    SV = [f]
    SV_deps = [()]
    AV = []
    FV = FVpre + [f]
    IM = {(): 0 + offset}
    compare_compute_argument_factorization(SV, SV_deps, AV, FV, IM)

    # Test basic non-argument sum
    SV = [f, g, f + g]
    SV_deps = [(), (), (0, 1)]
    AV = []
    FV = FVpre + [f, g, f + g]
    IM = {(): 2 + offset}
    compare_compute_argument_factorization(SV, SV_deps, AV, FV, IM)

    # Test basic non-argument product
    SV = [f, g, f * g]
    SV_deps = [(), (), (0, 1)]
    AV = []
    FV = FVpre + [f, g, f * g]
    IM = {(): 2 + offset}
    compare_compute_argument_factorization(SV, SV_deps, AV, FV, IM)

    # Test basic single-argument-only expression
    SV = [v]
    SV_deps = [()]
    AV = [v]
    FV = FVpre + []
    IM = {(0,): 1}  # v == AV[0] * FV[1]
    compare_compute_argument_factorization(SV, SV_deps, AV, FV, IM)

    # Test basic coefficient-argument product
    SV = [f, v, f * v]
    SV_deps = [(), (), (0, 1)]
    AV = [v]
    FV = FVpre + [f]
    IM = {(0,): offset}  # f*v == AV[0] * FV[1]
    compare_compute_argument_factorization(SV, SV_deps, AV, FV, IM)

    # Test basic argument product
    SV = [u, v, u * v]
    SV_deps = [(), (), (0, 1)]
    AV = [v, u]  # Test function < trial function
    FV = FVpre + []
    IM = {(0, 1): 1}  # v*u == (AV[0] * AV[1]) * FV[1]
    compare_compute_argument_factorization(SV, SV_deps, AV, FV, IM)

    # Test coefficient-argument products
    SV = [u, f, v, (f * v), u * (f * v)]
    SV_deps = [(), (), (), (1, 2), (0, 3)]
    AV = [v, u]
    FV = FVpre + [f]
    IM = {(0, 1): 0 + offset}  # f*(u*v) == (AV[0] * AV[1]) * FV[2]
    compare_compute_argument_factorization(SV, SV_deps, AV, FV, IM)

    # Test more complex situation
    SV = [u, u.dx(0), v,  # 0..2
          a, b, c, d, e,  # 3..7
          a * u, b * u.dx(0),  # 8..9
          c * v, d * v,  # 10..11
          a * u + b * u.dx(0),  # 12
          c * v + d * v,  # 13
          e * (a * u + b * u.dx(0)),  # 14
          (e * (a * u + b * u.dx(0))) * (c * v + d * v),  # 15
          ]
    SV_deps = [(), (), (),
                    (), (), (), (), (),
                    (0, 3), (1, 4),
                    (2, 5), (2, 6),
                    (8, 9),
                    (10, 11),
                    (7, 12),
                    (13, 14),
                    ]
    AV = [v, u, u.dx(0)]
    FV = FVpre + [a, b, c, d, e,  # 0..5
                  c + d,  # 6, introduced by SV[13]
                  e * a,  # 7, introduced by SV[14]
                  e * b,  # 8, introduced by SV[14]
                  (e * a) * (c + d),  # 9
                  (e * b) * (c + d),  # 10
                  ]
    IM = {(0, 1): 8 + offset,  # (a*e)*(c+d)*(u*v) == (AV[0] * AV[2]) * FV[13]
          (0, 2): 9 + offset}  # (b*e)*(c+d)*(u.dx(0)*v) == (AV[1] * AV[2]) * FV[12]
    compare_compute_argument_factorization(SV, SV_deps, AV, FV, IM)
