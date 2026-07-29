import basix.ufl
import numpy as np
import pytest
import ufl

import ffcx.codegeneration.jit
from ffcx.codegeneration.utils import dtype_to_c_type, dtype_to_scalar_dtype


@pytest.mark.parametrize("dtype", ["float64", "float32"])
@pytest.mark.parametrize("element", [("N1curl", {}), ("Lagrange", {"shape": (2,)})])
def test_interpolate(compile_args, dtype, element):

    cell = "triangle"
    family, el_kwargs = element

    domain = ufl.Mesh(basix.ufl.element("Lagrange", cell, 1, shape=(2,)))
    element = basix.ufl.element(family, cell, 2, **el_kwargs)
    V_int = ufl.FunctionSpace(domain, element)

    # Space containing other coefficients
    Q = ufl.FunctionSpace(domain, basix.ufl.element("Lagrange", cell, 2))
    w = ufl.Coefficient(Q)
    z = ufl.Coefficient(Q)
    q = ufl.Coefficient(Q)

    x = ufl.SpatialCoordinate(domain)
    c = ufl.Constant(domain)
    f = ufl.as_vector((z, -x[0])) + ufl.as_vector((q, q))
    If = ufl.Interpolate(f, V_int)
    J = (w * If[0] + If[1]) * ufl.dx
    f_ref = ufl.as_vector((x[1], -x[0])) + ufl.as_vector((1, 1))
    J_ref = (w * f_ref[0] + f_ref[1]) * ufl.dx

    compiled_forms, module, _code = ffcx.codegeneration.jit.compile_forms(
        [J, J_ref],
        options={"scalar_type": dtype},
        # cache_dir=".ffcx_cache",
        cffi_extra_compile_args=[],  # compile_args,
        visualise=False,
    )

    xdtype = dtype_to_scalar_dtype(dtype)
    scale = 4.2
    coords = np.array([[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, scale, 0.0]], dtype=xdtype).flatten()
    c = np.array([2.3], dtype=dtype)
    # Coefficients are ordered according to when they were created, thus
    # w, z, q
    # We set z to x[1], w, q to 1 and w to a some non-zero value
    q_size = Q.ufl_element().basix_element.dim
    assert q_size == 6
    d = np.empty(3 * q_size, dtype=dtype)
    d[0:q_size] = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7]  # w
    d[q_size : 2 * q_size] = [0.0, 0.0, scale * 1.0, 0.5 * scale, 0.5 * scale, 0.0]  # z
    d[2 * q_size : 3 * q_size] = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]  # q

    # Get kernel
    ffi = module.ffi
    form = compiled_forms[0]
    offsets = form.form_integral_offsets
    cell = module.lib.cell
    assert offsets[cell + 1] - offsets[cell] == 1

    default_integral = form.form_integrals[offsets[cell]]

    A = np.zeros(1, dtype=dtype)
    c_type, c_xtype = dtype_to_c_type(dtype), dtype_to_c_type(xdtype)
    kernel = getattr(default_integral, f"tabulate_tensor_{dtype}")
    kernel(
        ffi.cast(f"{c_type} *", A.ctypes.data),
        ffi.cast(f"{c_type} *", d.ctypes.data),
        ffi.cast(f"{c_type} *", c.ctypes.data),
        ffi.cast(f"{c_xtype} *", coords.ctypes.data),
        ffi.NULL,
        ffi.NULL,
        ffi.NULL,
    )

    d_ref = np.zeros(q_size, dtype=dtype)
    d_ref[:] = d[0:q_size]  # w
    c_ref = c.copy()
    form_ref = compiled_forms[1]
    offsets_ref = form_ref.form_integral_offsets
    assert offsets_ref[cell + 1] - offsets_ref[cell] == 1
    ref_integral = form_ref.form_integrals[offsets_ref[cell]]

    A_ref = np.zeros(1, dtype=dtype)
    ref_kernel = getattr(ref_integral, f"tabulate_tensor_{dtype}")
    ref_kernel(
        ffi.cast(f"{c_type} *", A_ref.ctypes.data),
        ffi.cast(f"{c_type} *", d_ref.ctypes.data),
        ffi.cast(f"{c_type} *", c_ref.ctypes.data),
        ffi.cast(f"{c_xtype} *", coords.ctypes.data),
        ffi.NULL,
        ffi.NULL,
        ffi.NULL,
    )
    tol = np.finfo(dtype).eps * 100
    np.testing.assert_allclose(A, A_ref, atol=tol)
