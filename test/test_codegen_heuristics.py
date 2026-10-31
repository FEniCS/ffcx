# Copyright (C) 2026 Garth N. Wells and Paul T. Kühner
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Tests for code-generation size and storage heuristics."""

import re

import basix.ufl
import pytest
import ufl

import ffcx.codegeneration.integral_generator as integral_generator
import ffcx.codegeneration.jit
import ffcx.codegeneration.lnodes as L
from ffcx.codegeneration.definitions import _should_unroll_coordinate_dofs
from ffcx.codegeneration.integral_generator import (
    _UNSUPPORTED_VECTOR_MATH_MIN_CONTRACTION_WORK,
    _contains_vectorizable_mathfunction,
    _contraction_work,
    _fw_cache_tile_size,
    _supports_vector_math_loop_split,
)
from ffcx.codegeneration.optimizer import (
    expr_key,
    prune_dead_scalars,
    reciprocal_cse_statements,
    scalar_declaration_names,
)


def test_expression_keys_distinguish_math_functions():
    x = L.Symbol("x", L.DataType.SCALAR)
    assert expr_key(L.MathFunction("sin", [x])) != expr_key(L.MathFunction("cos", [x]))


def test_coordinate_dof_unrolling_is_limited_to_small_non_tensor_geometry():
    """Curved simplex geometry must not be expanded into an unbounded sum."""
    assert _should_unroll_coordinate_dofs(10, False)
    assert not _should_unroll_coordinate_dofs(17, False)
    assert not _should_unroll_coordinate_dofs(4, True)


@pytest.mark.parametrize(
    "n_fw, num_points, expected",
    [
        (1, 10000, 2048),
        (256, 10000, 8),
        (257, 10000, 7),
        (2048, 10000, 1),
        (2049, 10000, None),
        (1000, 4, 2),
        (9, 11, 11),
        (81, 24, 24),
    ],
)
def test_fw_cache_tiles_do_not_exceed_the_stack_budget(n_fw, num_points, expected):
    """The cache must remain bounded even with many fw intermediates."""
    tile = _fw_cache_tile_size(n_fw, num_points)
    assert tile == expected
    if tile is not None:
        assert n_fw * tile * 16 <= 32768


def test_split_only_targets_supported_vector_math_hosts(monkeypatch):
    """The cache split is disabled where the compiler cannot use libmvec."""
    sin = L.MathFunction("sin", [L.Symbol("x", L.DataType.SCALAR)])
    power = L.MathFunction("power", [L.Symbol("x", L.DataType.SCALAR), L.LiteralInt(7)])
    assert _contains_vectorizable_mathfunction(sin)
    assert not _contains_vectorizable_mathfunction(power)

    monkeypatch.setattr(integral_generator.sys, "platform", "darwin")
    assert not _supports_vector_math_loop_split()
    monkeypatch.setattr(integral_generator.sys, "platform", "linux")
    monkeypatch.setattr(integral_generator.platform, "machine", lambda: "x86_64")
    assert _supports_vector_math_loop_split()


def test_contraction_work_distinguishes_small_and_large_contractions():
    """Unsupported targets split only when enough contraction work remains."""
    index = L.Symbol("i", L.DataType.INT)
    output = L.Symbol("output", L.DataType.SCALAR)
    small = L.ForRange(index, 0, 10, [L.Assign(output, 1.0) for _ in range(12)])
    large = L.ForRange(index, 0, 10, [L.Assign(output, 1.0) for _ in range(26)])
    assert _contraction_work(small) == 120
    assert _contraction_work(large) == 260
    assert _contraction_work(small) < _UNSUPPORTED_VECTOR_MATH_MIN_CONTRACTION_WORK
    assert _contraction_work(large) >= _UNSUPPORTED_VECTOR_MATH_MIN_CONTRACTION_WORK


def test_curved_simplex_geometry_uses_a_coordinate_dof_loop(compile_args):
    """A high-order simplex must not become one enormous coordinate sum."""
    coordinate_element = basix.ufl.element("Lagrange", "triangle", 5, shape=(2,))
    domain = ufl.Mesh(coordinate_element)
    element = basix.ufl.element("Lagrange", "triangle", 1)
    space = ufl.FunctionSpace(domain, element)
    test = ufl.TestFunction(space)
    x = ufl.SpatialCoordinate(domain)
    form = ufl.sin(x[0]) * test * ufl.dx(metadata={"quadrature_degree": 6})

    _, _, (_, implementation) = ffcx.codegeneration.jit.compile_forms(
        [form], options={"scalar_type": "float64"}, cffi_extra_compile_args=compile_args
    )
    assert "for (int ic = 0; ic < 21; ++ic)" in implementation


def test_tiled_factor_caches_keep_the_stack_budget(compile_args):
    """Factor caches of a tiled quadrature loop hold one tile, not all points."""
    element = basix.ufl.element("Lagrange", "tetrahedron", 2, shape=(3,))
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "tetrahedron", 1, shape=(3,)))
    space = ufl.FunctionSpace(domain, element)
    state = ufl.Coefficient(space)
    test, trial = ufl.TestFunction(space), ufl.TrialFunction(space)
    F = ufl.Identity(3) + ufl.grad(state)
    energy = ufl.inner(F, F) / 2 - ufl.ln(ufl.det(F)) + ufl.ln(ufl.det(F)) ** 2
    measure = ufl.dx(metadata={"quadrature_degree": 8})
    tangent = ufl.derivative(ufl.derivative(energy * measure, state, test), state, trial)

    _, _, (_, implementation) = ffcx.codegeneration.jit.compile_forms(
        [tangent], options={"scalar_type": "float64"}, cffi_extra_compile_args=compile_args
    )
    points = int(re.search(r"weights_\w+\[(\d+)\]", implementation)[1])
    sizes = [int(n) for n in re.findall(r"double fw\d+_q\[(\d+)\];", implementation)]
    assert sizes and all(size < points for size in sizes)
    assert len(sizes) * max(sizes) * 16 <= 32768


def test_dead_scalar_pruning_respects_section_declarations():
    """LNode liveness removes an unused declaration before C formatting."""
    live = L.Symbol("live", L.DataType.SCALAR)
    dead = L.Symbol("dead", L.DataType.SCALAR)
    output = L.Symbol("output", L.DataType.SCALAR)
    section = L.Section(
        "test",
        [L.Assign(output, live)],
        [L.VariableDecl(live, 1.0), L.VariableDecl(dead, 2.0)],
    )

    pruned = prune_dead_scalars([section], scalar_declaration_names([section]))
    assert [declaration.symbol.name for declaration in pruned[0].declarations] == ["live"]


@pytest.mark.parametrize("dtype", [L.DataType.REAL, L.DataType.SCALAR])
def test_reciprocal_cse(dtype):
    #            ==> r = 1 / c
    # x = a / c      x = a * r
    # y = b / c      y = b * r

    c = L.Symbol("c", dtype)
    exprs = [
        L.VariableDecl(L.Symbol("x", L.DataType.REAL), L.Div(L.Symbol("a", L.DataType.REAL), c)),
        L.VariableDecl(L.Symbol("y", dtype), L.Div(L.Symbol("b", dtype), c)),
    ]

    r, *exprs_recip = reciprocal_cse_statements(exprs)

    assert r.value == L.Div(L.LiteralFloat(1.0), c)
    assert all(expr.value.rhs.name == r.symbol.name for expr in exprs_recip)

    # (possibly) wider type of r dictates the expression dtype
    assert r.symbol.dtype == dtype
    assert all(statement.value.dtype == dtype for statement in exprs_recip)


def test_library_calls_leave_the_evaluation_loop(compile_args):
    """A library math call moves to a loop of its own, so evaluation can vectorise."""
    element = basix.ufl.element("Lagrange", "tetrahedron", 2, shape=(3,))
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "tetrahedron", 1, shape=(3,)))
    space = ufl.FunctionSpace(domain, element)
    state, test = ufl.Coefficient(space), ufl.TestFunction(space)
    F = ufl.Identity(3) + ufl.grad(state)
    energy = (ufl.inner(F, F) - 3) / 2 - ufl.ln(ufl.det(F)) + 0.625 * ufl.ln(ufl.det(F)) ** 2
    residual = ufl.derivative(energy * ufl.dx(metadata={"quadrature_degree": 4}), state, test)

    _, _, (_, implementation) = ffcx.codegeneration.jit.compile_forms(
        [residual], options={"scalar_type": "float64"}, cffi_extra_compile_args=compile_args
    )
    calls = re.findall(r"\blog\((\w+)", implementation)
    assert calls and all(argument.startswith("call_input") for argument in calls)
    loop = implementation.index("Evaluate and cache weighted scalar factors")
    assert "log(" not in implementation[loop:]


def test_outer_loops_keep_long_inner_loops_vectorised(compile_args):
    """Clang is told not to vectorise a quadrature loop across long inner loops."""
    element = basix.ufl.element("Lagrange", "tetrahedron", 3)
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "tetrahedron", 1, shape=(3,)))
    space = ufl.FunctionSpace(domain, element)
    f, v = ufl.Coefficient(space), ufl.TestFunction(space)
    short = ufl.FunctionSpace(domain, basix.ufl.element("Lagrange", "tetrahedron", 1))
    g, w = ufl.Coefficient(short), ufl.TestFunction(short)

    def source(form):
        _, _, (_, implementation) = ffcx.codegeneration.jit.compile_forms(
            [form], options={"scalar_type": "float64"}, cffi_extra_compile_args=compile_args
        )
        body = implementation[implementation.index("void tabulate_tensor") :]
        return re.sub(r"\n\s*", "\n", body)

    # Twenty dofs: the point loop stays scalar so the dof loops vectorise
    assert "FFCX_NO_VECTORIZE\nfor (int iq" in source(f * v * ufl.dx)
    # Four dofs: vectorising across points remains Clang's choice
    assert "FFCX_NO_VECTORIZE" not in source(g * w * ufl.dx)
