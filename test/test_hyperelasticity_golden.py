# Copyright (C) 2026 Garth N. Wells
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Golden-master numeric regression test for demo/HyperElasticity.py.

This form (P2 vector displacement, nonlinear strain energy, Jacobian built via
``ufl.derivative()``) is the fixture an investigation into FFCx's factorization-graph
size for deeply-nested tensor composition is anchored to. It previously had no numeric
correctness coverage at all -- ``demo/test_demos.py`` only checks that the demo
compiles, never asserts a result -- so any future change to how nested tensor
expressions get scalarized/factorized had nothing to catch a silent wrong-value
regression. This test freezes the current compiler's output for both the residual
(``a_F``) and Jacobian (``a_J``) cell integrals against fixed, deterministic,
non-degenerate input data, so any future change to this pipeline can be checked against
these values.

The frozen data was generated once with the unmodified compiler; see
``test/data/hyperelasticity_golden.npz`` (coords, w, c, a_F, a_J).
"""

import importlib.util
import types
from pathlib import Path

import numpy as np
import pytest

import ffcx.codegeneration.jit

demo_dir = Path(__file__).parent.parent / "demo"
golden_file = Path(__file__).parent / "data" / "hyperelasticity_golden.npz"


def _load_hyperelasticity_demo() -> types.ModuleType:
    """Load demo/HyperElasticity.py as a module without touching sys.path."""
    spec = importlib.util.spec_from_file_location(
        "HyperElasticity", demo_dir / "HyperElasticity.py"
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize(
    "scalar_type, scalar_dtype, geometry_dtype, c_scalar_type, c_geometry_type, rtol, atol",
    [
        ("float32", np.float32, np.float32, "float", "float", 2e-5, 2e-6),
        ("float64", np.float64, np.float64, "double", "double", 1e-12, 1e-14),
        ("complex64", np.complex64, np.float32, "float _Complex", "float", 2e-5, 2e-6),
        (
            "complex128",
            np.complex128,
            np.float64,
            "double _Complex",
            "double",
            1e-12,
            1e-14,
        ),
    ],
)
def test_hyperelasticity_cell_kernels_match_golden_values(
    compile_args,
    scalar_type,
    scalar_dtype,
    geometry_dtype,
    c_scalar_type,
    c_geometry_type,
    rtol,
    atol,
) -> None:
    """Tabulate demo/HyperElasticity.py's a_F and a_J cell kernels, compare to frozen values."""
    demo = _load_hyperelasticity_demo()
    golden = np.load(golden_file)

    compiled_forms, module, _code = ffcx.codegeneration.jit.compile_forms(
        demo.forms, options={"scalar_type": scalar_type}, cffi_extra_compile_args=compile_args
    )
    ffi = module.ffi
    cell = module.lib.cell

    coords = golden["coords"].astype(geometry_dtype)
    w = golden["w"].astype(scalar_dtype)
    c = golden["c"].astype(scalar_dtype)

    for name, compiled_form, size in [
        ("a_F", compiled_forms[0], 30),
        ("a_J", compiled_forms[1], 900),
    ]:
        offsets = compiled_form.form_integral_offsets
        integral = compiled_form.form_integrals[offsets[cell]]
        A = np.zeros(size, dtype=scalar_dtype)
        kernel = getattr(integral, f"tabulate_tensor_{scalar_type}")
        kernel(
            ffi.cast(f"{c_scalar_type} *", A.ctypes.data),
            ffi.cast(f"{c_scalar_type} *", w.ctypes.data),
            ffi.cast(f"{c_scalar_type} *", c.ctypes.data),
            ffi.cast(f"{c_geometry_type} *", coords.ctypes.data),
            ffi.NULL,
            ffi.NULL,
            ffi.NULL,
        )
        np.testing.assert_allclose(
            A,
            golden[name].astype(scalar_dtype),
            rtol=rtol,
            atol=atol,
            err_msg=f"HyperElasticity {name} cell kernel no longer matches the frozen golden "
            "value -- if this change is intentional (e.g. a deliberate scalarization/"
            "factorization rewrite), regenerate test/data/hyperelasticity_golden.npz and "
            "confirm the new values via an independent method before updating it.",
        )
