# Copyright (C) 2025 Paul T. Kühner
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import ctypes
import importlib.util
from pathlib import Path
from types import ModuleType

import numpy as np
import numpy.typing as npt
import pytest

import ffcx.main
from ffcx.codegeneration.numba.utils import ufcx_kernel_signature
from ffcx.codegeneration.utils import dtype_to_scalar_dtype

numba = pytest.importorskip("numba")


def wrap_kernel(scalar_type, real_type):
    c_signature = ufcx_kernel_signature(scalar_type, real_type)
    return numba.cfunc(c_signature, nopython=True)


def as_C_array(np_array: npt.NDArray):
    dtype_C = np.ctypeslib.as_ctypes_type(np_array.dtype)
    return np_array.ctypes.data_as(ctypes.POINTER(dtype_C))


def load_module(path: Path) -> ModuleType:
    """Load a generated module from an explicit file path."""
    spec = importlib.util.spec_from_file_location(path.stem, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_external_module(tmp_path: Path):
    """Inject ffcx.codegeneration.numba as external module, to test support for custom backends."""
    opts = "--language ffcx.codegeneration.numba --scalar_type float64"
    dir = Path(__file__).parent
    args = [str(dir / "poisson.py"), *opts.split(" "), "--dir", str(tmp_path)]
    assert ffcx.main.main(args) == 0


@pytest.fixture(scope="module")
def generate_poisson(tmp_path_factory: pytest.TempPathFactory) -> ModuleType:
    """Generate the numba Poisson module and load it.

    The module is generated into a pytest temporary directory and loaded
    from that path, rather than written to and imported by name from the
    current working directory. A shared file in the working directory is
    racy when tests run in parallel (``pytest-xdist``): one worker can
    import the module while another is still writing it.

    Module-scoped so the generation runs once rather than once per test.
    The tests only read from the module, so sharing one instance is safe.
    ``tmp_path_factory`` is used because the function-scoped ``tmp_path``
    cannot be requested from a module-scoped fixture.
    """
    tmp_path = tmp_path_factory.mktemp("numba")
    dir = Path(__file__).parent
    args = [str(dir / "poisson.py"), "--language", "numba", "--dir", str(tmp_path)]
    assert ffcx.main.main(args) == 0
    return load_module(tmp_path / "poisson_numba.py")


@pytest.mark.parametrize("scalar_type", ["float32", "float64"])  # TODO: complex limited by ctypes
def test_integral(scalar_type: str, generate_poisson: ModuleType) -> None:
    poisson = generate_poisson

    dtype = np.dtype(scalar_type).type
    dtype_r = dtype_to_scalar_dtype(dtype)

    kernel_a = wrap_kernel(dtype, dtype_r)(poisson.form_poisson_a.form_integrals[0].tabulate_tensor)

    A = np.zeros((3, 3), dtype=dtype)
    w = np.array([], dtype=dtype)
    kappa_value = np.array([[1.0, 2.0], [3.0, 4.0]])
    c = np.array(kappa_value.flatten(), dtype=dtype)
    coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=dtype_r)
    empty = np.empty((0,), dtype=dtype_r)

    kernel_a(
        as_C_array(A),
        as_C_array(w),
        as_C_array(c),
        as_C_array(coords),
        as_C_array(empty),
        as_C_array(empty),
        0,
    )
    A_expected = np.array(
        [[1.0, -0.5, -0.5], [-0.5, 0.5, 0.0], [-0.5, 0.0, 0.5]], dtype=scalar_type
    )

    assert np.allclose(A, np.trace(kappa_value) * A_expected)

    kernel_L = wrap_kernel(dtype, dtype_r)(poisson.form_poisson_L.form_integrals[0].tabulate_tensor)

    b = np.zeros((3,), dtype=dtype)
    w = np.full((3,), 0.5, dtype=dtype)
    c = np.empty((0,), dtype=dtype)
    coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=dtype_r)
    empty = np.empty((0,), dtype=dtype_r)

    kernel_L(
        as_C_array(b),
        as_C_array(w),
        as_C_array(c),
        as_C_array(coords),
        as_C_array(empty),
        as_C_array(empty),
        0,
    )

    b_expected = np.full((3,), 1 / 6, dtype=np.float64)
    assert np.allclose(b, 0.5 * b_expected)


@pytest.mark.parametrize("scalar_type", ["float32", "float64"])  # TODO: complex limited by ctypes
def test_expression(scalar_type: str, generate_poisson: ModuleType) -> None:
    poisson = generate_poisson

    dtype = np.dtype(scalar_type).type
    dtype_r = dtype_to_scalar_dtype(dtype)

    kernel_expr = wrap_kernel(dtype, dtype_r)(poisson.expression_poisson_0.tabulate_tensor)

    e = np.zeros((6 * 3,), dtype=dtype)
    w = np.array(
        [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11], dtype=dtype
    )
    kappa_value = np.array([[1.0, 2.0], [3.0, 4.0]])
    c = np.array(kappa_value.flatten(), dtype=dtype)
    coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=dtype_r)
    empty = np.empty((0,), dtype=dtype_r)

    kernel_expr(
        as_C_array(e),
        as_C_array(w),
        as_C_array(c),
        as_C_array(coords),
        as_C_array(empty),
        as_C_array(empty),
        0,
    )
    e_expected = np.array(
        [5, 7, 8, 11, 15, 18, 14, 16, 17, 32, 36, 39, 23, 25, 26, 53, 57, 60], dtype=dtype
    )
    assert np.allclose(e, e_expected)
