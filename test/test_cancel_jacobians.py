import time

import basix.ufl
import numpy as np
import pytest
import ufl

import ffcx.codegeneration.jit
from ffcx.codegeneration.utils import dtype_to_c_type, dtype_to_scalar_dtype


@pytest.mark.parametrize(
    "dtype",
    [
        "float32",
        "float64",
    ],
)
@pytest.mark.parametrize("degree", [1, 3, 5])
@pytest.mark.parametrize(
    "cell_type, gdim",
    [("triangle", 2), ("tetrahedron", 3), ("quadrilateral", 2), ("hexahedron", 3)],
)
def test_cancel_jacobians(dtype, compile_args, degree, cell_type, gdim):
    basix.cell
    domain = ufl.Mesh(basix.ufl.element("Lagrange", cell_type, 1, shape=(gdim,)))
    rt = basix.ufl.element("Raviart-Thomas", cell_type, degree)
    space = ufl.FunctionSpace(domain, rt)
    u, v = ufl.TrialFunction(space), ufl.TestFunction(space)
    form = ufl.inner(ufl.div(u), ufl.div(v)) * ufl.dx

    forms = [form]
    compiled_forms, module, _code = ffcx.codegeneration.jit.compile_forms(
        forms,
        options={"scalar_type": dtype, "do_cancel_jacobian_products": True},
        cffi_extra_compile_args=compile_args,
    )

    A = np.zeros((rt.dim, rt.dim), dtype=dtype)
    w = np.array([], dtype=dtype)
    c = np.array([], dtype=dtype)

    xdtype = dtype_to_scalar_dtype(dtype)
    if cell_type == "triangle":
        coords = np.array([[0.0, 0.0, 0.0], [2.0, 0.1, 0.0], [0.3, 1.1, 0.0]], dtype=xdtype)
    elif cell_type == "quadrilateral":
        coords = np.array(
            [[0.0, 0.0, 0.0], [1.1, 0.0, 0.0], [-0.1, 1.3, 0.0], [1.2, 1.2, 0.0]], dtype=xdtype
        )
    elif cell_type == "tetrahedron":
        coords = np.array(
            [[0.0, 0.0, 0.0], [1.1, 0.0, 0.0], [0.0, 1.2, 0.0], [-0.1, -0.2, 0.8]], dtype=xdtype
        )
    elif cell_type == "hexahedron":
        coords = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.1, 0.0, 0.0],
                [0.0, 0.9, 0.0],
                [1.2, 1.1, 0.0],
                [0.0, 0.0, 0.8],
                [1.0, 0.0, 1.0],
                [0.0, 1.0, 1.0],
                [1.0, 1.0, 1.0],
            ],
            dtype=xdtype,
        )
    else:
        raise ValueError(f"Unsupported {cell_type=}")

    c_type, c_xtype = dtype_to_c_type(dtype), dtype_to_c_type(xdtype)
    ffi = module.ffi
    form = compiled_forms[0]
    default_integral = form.form_integrals[0]
    kernel = getattr(default_integral, f"tabulate_tensor_{dtype}")

    new_start = time.perf_counter()
    kernel(
        ffi.cast(f"{c_type} *", A.ctypes.data),
        ffi.cast(f"{c_type} *", w.ctypes.data),
        ffi.cast(f"{c_type} *", c.ctypes.data),
        ffi.cast(f"{c_xtype} *", coords.ctypes.data),
        ffi.NULL,
        ffi.NULL,
        ffi.NULL,
    )
    new_end = time.perf_counter()

    old_compiled_forms, old_module, _code_old = ffcx.codegeneration.jit.compile_forms(
        forms,
        options={"scalar_type": dtype, "do_cancel_jacobian_products": False},
        cffi_extra_compile_args=compile_args,
    )

    A_ref = np.zeros((rt.dim, rt.dim), dtype=dtype)

    c_type, c_xtype = dtype_to_c_type(dtype), dtype_to_c_type(xdtype)
    ffi = old_module.ffi
    form = old_compiled_forms[0]
    default_integral = form.form_integrals[0]
    kernel = getattr(default_integral, f"tabulate_tensor_{dtype}")

    old_start = time.perf_counter()
    kernel(
        ffi.cast(f"{c_type} *", A_ref.ctypes.data),
        ffi.cast(f"{c_type} *", w.ctypes.data),
        ffi.cast(f"{c_type} *", c.ctypes.data),
        ffi.cast(f"{c_xtype} *", coords.ctypes.data),
        ffi.NULL,
        ffi.NULL,
        ffi.NULL,
    )
    old_end = time.perf_counter()

    np.testing.assert_allclose(
        A, A_ref, rtol=np.finfo(dtype).eps * 1000, atol=np.finfo(dtype).eps * 1000
    )

    print(f"Time with cancel_jacobian_products: {new_end - new_start:.6e} seconds")
    print(f"Time without cancel_jacobian_products: {old_end - old_start:.6e} seconds")
    print(f"Speedup: {(old_end - old_start) / (new_end - new_start):.5f}x")
