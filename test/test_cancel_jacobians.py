import basix.ufl
import numpy as np
import pytest
import ufl

import ffcx.codegeneration.jit
from ffcx.codegeneration.utils import dtype_to_c_type, dtype_to_scalar_dtype


def _tabulate_tensor(forms, shape, dtype, coords, compile_args, do_cancel_jacobian_products=None):
    """Compile and tabulate the default integral of a form."""
    options = {"scalar_type": dtype}
    if do_cancel_jacobian_products is not None:
        options["do_cancel_jacobian_products"] = do_cancel_jacobian_products
    compiled_forms, module, _ = ffcx.codegeneration.jit.compile_forms(
        forms,
        options=options,
        cffi_extra_compile_args=compile_args,
    )

    tensor = np.zeros(shape, dtype=dtype)
    w = np.array([], dtype=dtype)
    c = np.array([], dtype=dtype)
    xdtype = dtype_to_scalar_dtype(dtype)
    c_type, c_xtype = dtype_to_c_type(dtype), dtype_to_c_type(xdtype)

    default_integral = compiled_forms[0].form_integrals[0]
    kernel = getattr(default_integral, f"tabulate_tensor_{dtype}")
    ffi = module.ffi
    kernel(
        ffi.cast(f"{c_type} *", tensor.ctypes.data),
        ffi.cast(f"{c_type} *", w.ctypes.data),
        ffi.cast(f"{c_type} *", c.ctypes.data),
        ffi.cast(f"{c_xtype} *", coords.ctypes.data),
        ffi.NULL,
        ffi.NULL,
        ffi.NULL,
    )

    return tensor


def _triangle_coordinates(dtype):
    """Return non-degenerate triangle coordinates in the scalar real dtype."""
    xdtype = dtype_to_scalar_dtype(dtype)
    return np.array([[0.0, 0.0, 0.0], [2.0, 0.1, 0.0], [0.3, 1.1, 0.0]], dtype=xdtype)


def _assert_tensors_equal(tensor, reference, dtype):
    """Check tensor equivalence with a tolerance appropriate for the scalar type."""
    if dtype in ("float32", "complex64"):
        # The algebraically equivalent paths evaluate J @ inv(J) in a different
        # order, so the uncancelled path accumulates more single-precision round-off.
        tol = 5e-4
    elif dtype in ("float64", "complex128"):
        tol = 1e-12
    else:
        raise ValueError(f"Unsupported {dtype=}")
    np.testing.assert_allclose(tensor, reference, rtol=tol, atol=tol)


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
    [
        ("triangle", 2),
        ("triangle", 3),
        ("tetrahedron", 3),
        ("quadrilateral", 2),
        ("hexahedron", 3),
    ],
)
def test_cancel_jacobians(dtype, compile_args, degree, cell_type, gdim):
    domain = ufl.Mesh(basix.ufl.element("Lagrange", cell_type, 1, shape=(gdim,)))
    rt = basix.ufl.element("Raviart-Thomas", cell_type, degree)
    space = ufl.FunctionSpace(domain, rt)

    u, v = ufl.TrialFunction(space), ufl.TestFunction(space)
    form = ufl.inner(ufl.div(u), ufl.div(v)) * ufl.dx

    forms = [form]
    xdtype = dtype_to_scalar_dtype(dtype)
    if cell_type == "triangle":
        coords = _triangle_coordinates(dtype)
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

    A = _tabulate_tensor(forms, (rt.dim, rt.dim), dtype, coords, compile_args)
    A_ref = _tabulate_tensor(forms, (rt.dim, rt.dim), dtype, coords, compile_args, False)
    _assert_tensors_equal(A, A_ref, dtype)


@pytest.mark.parametrize("dtype", ["complex64", "complex128"])
def test_cancel_jacobians_complex(dtype, compile_args):
    """Test cancellation for complex-valued RT div-div kernels."""
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
    rt = basix.ufl.element("Raviart-Thomas", "triangle", 3)
    space = ufl.FunctionSpace(domain, rt)
    u, v = ufl.TrialFunction(space), ufl.TestFunction(space)
    form = ufl.inner(ufl.div(u), ufl.div(v)) * ufl.dx

    coords = _triangle_coordinates(dtype)
    A = _tabulate_tensor([form], (rt.dim, rt.dim), dtype, coords, compile_args)
    A_ref = _tabulate_tensor([form], (rt.dim, rt.dim), dtype, coords, compile_args, False)
    _assert_tensors_equal(A, A_ref, dtype)


@pytest.mark.parametrize("dtype", ["float64", "complex128"])
def test_cancel_both_jacobian_product_orders(dtype, compile_args):
    """Test both J @ inv(J) and inv(J) @ J contractions."""
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
    element = basix.ufl.element("Lagrange", "triangle", 1)
    space = ufl.FunctionSpace(domain, element)
    u, v = ufl.TrialFunction(space), ufl.TestFunction(space)
    J = ufl.Jacobian(domain)
    K = ufl.JacobianInverse(domain)
    i, k = ufl.indices(2)
    jacobian_products = J[i, k] * K[k, i] + K[i, k] * J[k, i]
    form = jacobian_products * ufl.inner(u, v) * ufl.dx

    coords = _triangle_coordinates(dtype)
    A = _tabulate_tensor([form], (element.dim, element.dim), dtype, coords, compile_args)
    A_ref = _tabulate_tensor([form], (element.dim, element.dim), dtype, coords, compile_args, False)
    _assert_tensors_equal(A, A_ref, dtype)
