# Copyright (C) 2024 Igor A. Baratta
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import sys

import basix.ufl
import cffi
import numpy as np
import pytest
import ufl

import ffcx.codegeneration.jit
import ffcx.codegeneration.utils as utils


def generate_kernel(forms, scalar_type, options):
    """Generate kernel for given forms."""
    compiled_forms, module, code = ffcx.codegeneration.jit.compile_forms(
        forms, options={"scalar_type": scalar_type}
    )

    for f, compiled_f in zip(forms, compiled_forms):
        assert compiled_f.rank == len(f.arguments())

    form0 = compiled_forms[0]

    offsets = form0.form_integral_offsets
    cell = module.lib.cell
    assert offsets[cell + 1] - offsets[cell] == 1
    integral_id = form0.form_integral_ids[offsets[cell]]
    assert integral_id == -1
    default_integral = form0.form_integrals[offsets[cell]]
    kernel = getattr(default_integral, f"tabulate_tensor_{scalar_type}")
    return kernel, code, module


@pytest.mark.parametrize(
    "dtype",
    [
        "float32",
        "float64",
        pytest.param(
            "complex64",
            marks=pytest.mark.xfail(
                sys.platform.startswith("win32"),
                raises=NotImplementedError,
                reason="missing _Complex",
            ),
        ),
        pytest.param(
            "complex128",
            marks=pytest.mark.xfail(
                sys.platform.startswith("win32"),
                raises=NotImplementedError,
                reason="missing _Complex",
            ),
        ),
    ],
)
def test_numba_kernel_signature(dtype):
    try:
        import numba
    except ImportError:
        pytest.skip("Numba not installed")

    # Create a simple form
    mesh = ufl.Mesh(basix.ufl.element("P", "triangle", 2, shape=(2,)))
    e = basix.ufl.element("Lagrange", "triangle", 2)

    V = ufl.FunctionSpace(mesh, e)
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)

    a = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx

    # Generate and compile the kernel
    kernel, _code, _module = generate_kernel([a], dtype, {})

    # Convert to numpy dtype
    np_dtype = np.dtype(dtype)

    # Generate the Numba signature
    xtype = utils.dtype_to_scalar_dtype(dtype)
    signature = utils.numba_ufcx_kernel_signature(np_dtype, xtype)
    assert isinstance(signature, numba.core.typing.templates.Signature)

    # Get the signature from the compiled kernel
    ffi = cffi.FFI()
    args = ffi.typeof(kernel).args

    # check that the signature is equivalent to the one in the generated code
    assert len(args) == len(signature.args)
    for i, (arg, sig) in enumerate(zip(args, signature.args)):
        type_name = sig.name.replace(str(np_dtype), utils.dtype_to_c_type(np_dtype))
        ctypes_name = type_name.replace(" *", "*")
        assert ctypes_name == type_name
