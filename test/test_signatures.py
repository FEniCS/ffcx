# Copyright (C) 2024 Igor A. Baratta
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import basix.ufl
import cffi
import numpy as np
import pytest
import ufl

import ffcx.codegeneration.jit
import ffcx.codegeneration.utils as utils


def generate_kernel(forms, scalar_type, options):
    """Generate kernel for given forms."""
    compiled_forms, module, code = ffcx.codegeneration.jit.compile_forms(forms)

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


@pytest.mark.parametrize("dtype", [np.float32, np.float64, np.complex64, np.complex128])
def test_numba_kernel_signature(dtype):
    try:
        import numba
    except ImportError:
        pytest.skip("Numba not installed")

    # Convert to numpy dtype
    dtype = np.dtype(dtype)

    # Create a simple form
    mesh = ufl.Mesh(basix.ufl.element("P", "triangle", 2, shape=(2,)))
    e = basix.ufl.element("Lagrange", "triangle", 2)

    V = ufl.FunctionSpace(mesh, e)
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)

    a = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx

    # Generate and compile the kernel
    kernel, code, module = generate_kernel([a], dtype, {})

    # Generate the Numba signature
    xtype = utils.dtype_to_scalar_dtype(dtype)
    signature = utils.numba_ufcx_kernel_signature(dtype, xtype)
    assert isinstance(signature, numba.core.typing.templates.Signature)

    # Get the signature from the compiled kernel
    ffi = cffi.FFI()
    args = ffi.typeof(kernel).args

    # check that the signature is equivalent to the one in the generated code
    assert len(args) == len(signature.args)
    for i, (arg, sig) in enumerate(zip(args, signature.args)):
        type_name = sig.name.replace(str(dtype), utils.dtype_to_c_type(dtype))
        ctypes_name = type_name.replace(" *", "*")
        assert ctypes_name == type_name
