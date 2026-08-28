# Copyright (C) 2026 Garth N. Wells
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Benchmark suite of forms for measuring generated cell-kernel performance.

Each :class:`BenchCase` pairs a UFL form with sample input data (coordinates,
coefficients, constants) for its single cell integral, so ``bench/harness.py``
can compile it and time repeated calls to the generated ``tabulate_tensor``.
"""

import importlib.util
from pathlib import Path

import basix.ufl
import numpy as np
import ufl

_HERE = Path(__file__).parent
_DEMO_DIR = _HERE.parent / "demo"
_GOLDEN_FILE = _HERE.parent / "test" / "data" / "hyperelasticity_golden.npz"

_REFERENCE_VERTICES = {
    "triangle": [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)],
    "tetrahedron": [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)],
}


def _reference_coords(cellname: str) -> np.ndarray:
    """Flattened vertex coordinates of the reference cell (P1 geometry)."""
    return np.array(_REFERENCE_VERTICES[cellname], dtype=np.float64).flatten()


class BenchCase:
    """A single (form, sample input data) case for the benchmark harness."""

    def __init__(
        self,
        name: str,
        form: ufl.Form,
        cellname: str,
        coords: np.ndarray,
        w: np.ndarray | None = None,
        c: np.ndarray | None = None,
    ):
        """Initialise a benchmark case."""
        self.name = name
        self.form = form
        self.cellname = cellname
        self.coords = coords
        self.w = w if w is not None else np.zeros(0, dtype=np.float64)
        self.c = c if c is not None else np.zeros(0, dtype=np.float64)


def _scalar_cases(cellname: str, gdim: int, degree: int) -> list[BenchCase]:
    """Mass and Poisson forms for a scalar Lagrange element."""
    element = basix.ufl.element("Lagrange", cellname, degree)
    domain = ufl.Mesh(basix.ufl.element("Lagrange", cellname, 1, shape=(gdim,)))
    space = ufl.FunctionSpace(domain, element)
    u, v = ufl.TrialFunction(space), ufl.TestFunction(space)
    coords = _reference_coords(cellname)
    return [
        BenchCase(f"mass_P{degree}_{cellname}", u * v * ufl.dx, cellname, coords),
        BenchCase(
            f"poisson_P{degree}_{cellname}",
            ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx,
            cellname,
            coords,
        ),
    ]


def _elasticity_cases(cellname: str, gdim: int, degree: int) -> list[BenchCase]:
    """Linear-elasticity bilinear form for a vector Lagrange element."""
    element = basix.ufl.element("Lagrange", cellname, degree, shape=(gdim,))
    domain = ufl.Mesh(basix.ufl.element("Lagrange", cellname, 1, shape=(gdim,)))
    space = ufl.FunctionSpace(domain, element)
    u, v = ufl.TrialFunction(space), ufl.TestFunction(space)
    coords = _reference_coords(cellname)
    form = ufl.inner(ufl.sym(ufl.grad(u)), ufl.sym(ufl.grad(v))) * ufl.dx
    return [BenchCase(f"elasticity_P{degree}_{cellname}", form, cellname, coords)]


def _pin_cell_quadrature_degree(form: ufl.Form, degree: int) -> ufl.Form:
    """Return a copy of ``form`` with every cell integral's quadrature degree pinned."""
    integrals = [
        i.reconstruct(metadata={"quadrature_degree": degree}) if i.integral_type() == "cell" else i
        for i in form.integrals()
    ]
    return ufl.Form(integrals)


def _hyperelasticity_cases() -> list[BenchCase]:
    """Residual and Jacobian cell kernels from demo/HyperElasticity.py.

    Reuses the frozen coordinates/coefficients/constants from the golden-master
    regression test (test/data/hyperelasticity_golden.npz) so the sample data is
    deterministic and non-degenerate.
    """
    spec = importlib.util.spec_from_file_location(
        "HyperElasticity", _DEMO_DIR / "HyperElasticity.py"
    )
    assert spec is not None and spec.loader is not None
    demo = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(demo)

    golden = np.load(_GOLDEN_FILE)
    coords, w, c = golden["coords"], golden["w"], golden["c"]

    a_F, a_J = demo.forms
    cases = []
    for name, form in [("hyperelasticity_a_F", a_F), ("hyperelasticity_a_J", a_J)]:
        cases.append(BenchCase(f"{name}_default_degree", form, "tetrahedron", coords, w, c))
        # The default UFL degree estimate for this form is very high (2744-point
        # rule on a tetrahedron -- see the perf-stack plan); a pinned low degree
        # is included so quadrature-count cost can be separated from per-point
        # kernel cost.
        cases.append(
            BenchCase(
                f"{name}_pinned_degree4",
                _pin_cell_quadrature_degree(form, 4),
                "tetrahedron",
                coords,
                w,
                c,
            )
        )
    return cases


def suite() -> list[BenchCase]:
    """Full benchmark suite."""
    cases: list[BenchCase] = []
    for cellname, gdim in [("triangle", 2), ("tetrahedron", 3)]:
        for degree in (1, 2):
            cases += _scalar_cases(cellname, gdim, degree)
        cases += _elasticity_cases(cellname, gdim, 1)
        cases += _elasticity_cases(cellname, gdim, 2)
    cases += _hyperelasticity_cases()
    return cases
