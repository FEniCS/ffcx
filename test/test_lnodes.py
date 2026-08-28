import importlib

import basix.ufl
import numpy as np
import pytest
import ufl
from cffi import FFI

import ffcx.codegeneration.jit
from ffcx.codegeneration import lnodes as L
from ffcx.codegeneration.C import Formatter
from ffcx.codegeneration.utils import dtype_to_c_type


def _power_lnode(exponent):
    """Translate ``base**exponent`` through ufl_to_lnodes, given a symbolic base.

    Uses ufl.variable to stop UFL constant-folding the power away before it
    reaches ufl_to_lnodes.
    """
    x = ufl.variable(ufl.constantvalue.FloatValue(3.0))
    op = x**exponent
    base = L.Symbol("x", L.DataType.SCALAR)
    exponent_lnode = L.ufl_to_lnodes(op.ufl_operands[1])
    return L.ufl_to_lnodes(op, base, exponent_lnode)


@pytest.mark.parametrize(
    "exponent,expected_type",
    # Exponent 0 and 1 aren't tested: UFL itself constant-folds/simplifies
    # x**0 and x**1 away before a Power node ever reaches ufl_to_lnodes.
    [(2, L.Mul), (3, L.Mul), (-1, L.Div), (-2, L.Div)],
)
def test_power_expansion_small_integer_exponents(exponent, expected_type):
    """Small integer powers expand to multiplications/a reciprocal, not pow()."""
    result = _power_lnode(exponent)
    assert isinstance(result, expected_type)
    assert not isinstance(result, L.MathFunction)


@pytest.mark.parametrize("exponent", [5, -5, 0.5])
def test_power_expansion_falls_back_to_pow_outside_small_range(exponent):
    """Exponents outside the small integer range still emit a pow() call."""
    result = _power_lnode(exponent)
    assert isinstance(result, L.MathFunction)
    assert result.function == "power"


def test_squared_coefficient_generates_no_pow_call(compile_args):
    """A squared term in a form compiles without a pow() call in the generated code.

    Also matters beyond the arithmetic: a stray pow() call makes
    integral_generator._contains_mathfunction treat the kernel as containing
    a transcendental call, triggering the quadrature-loop SIMD split for a
    form that has no actual transcendental function to vectorise.
    """
    element = basix.ufl.element("Lagrange", "triangle", 2)
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
    space = ufl.FunctionSpace(domain, element)
    u = ufl.Coefficient(space)
    v = ufl.TestFunction(space)
    form = (u**2) * v * ufl.dx

    _, _, (_, impl) = ffcx.codegeneration.jit.compile_forms(
        [form], options={"scalar_type": "float64"}, cffi_extra_compile_args=compile_args
    )
    assert "pow(" not in impl


@pytest.mark.parametrize("dtype", ("float32", "float64", "intc"))
def test_gemm(dtype):
    # Test LNodes simple matrix-matrix multiply in C
    p, q, r = 5, 16, 12

    A = L.Symbol("A", dtype=L.DataType.SCALAR)
    B = L.Symbol("B", dtype=L.DataType.SCALAR)
    C = L.Symbol("C", dtype=L.DataType.SCALAR)
    code = [L.Comment(f"Matrix multiply A{p, r} = B{p, q} * C{q, r}")]

    i = L.Symbol("i", dtype=L.DataType.INT)
    j = L.Symbol("j", dtype=L.DataType.INT)
    k = L.Symbol("k", dtype=L.DataType.INT)
    m_ij = L.MultiIndex([i, j], [p, q])
    m_ik = L.MultiIndex([i, k], [p, r])
    m_jk = L.MultiIndex([j, k], [q, r])

    body = [L.AssignAdd(A[m_ik], B[m_ij] * C[m_jk])]
    body = [L.ForRange(i, 0, p, body=body)]
    body = [L.ForRange(j, 0, q, body=body)]
    code += [L.ForRange(k, 0, r, body=body)]

    # Format into C and compile with CFFI
    format = Formatter(dtype=dtype)
    c_scalar = dtype_to_c_type(dtype)
    decl = f"void gemm({c_scalar} *A, {c_scalar} *B, {c_scalar} *C)"
    c_code = decl + "{\n" + format(L.StatementList(code)) + "\n}\n"

    ffibuilder = FFI()
    ffibuilder.cdef(decl + ";")
    ffibuilder.set_source(f"_gemm_{c_scalar}", c_code)
    ffibuilder.compile(verbose=True)
    _gemm = importlib.import_module(f"_gemm_{c_scalar}")
    gemm = _gemm.lib.gemm
    ffi = _gemm.ffi

    A = np.zeros((p, r), dtype=dtype)
    B = np.ones((p, q), dtype=dtype)
    C = np.ones((q, r), dtype=dtype)
    pA = ffi.cast(f"{c_scalar} *", A.ctypes.data)
    pB = ffi.cast(f"{c_scalar} *", B.ctypes.data)
    pC = ffi.cast(f"{c_scalar} *", C.ctypes.data)

    gemm(pA, pB, pC)
    assert np.all(A == q)


@pytest.mark.parametrize("dtype", ("float32", "float64", "intc"))
def test_gemv(dtype):
    # Test LNodes simple matvec multiply in C
    p, q = 5, 16

    y = L.Symbol("y", dtype=L.DataType.SCALAR)
    A = L.Symbol("A", dtype=L.DataType.SCALAR)
    x = L.Symbol("x", dtype=L.DataType.SCALAR)
    code = [L.Comment(f"Matrix-vector multiply y({p}) = A{p, q} * x({q})")]

    i = L.Symbol("i", dtype=L.DataType.INT)
    j = L.Symbol("j", dtype=L.DataType.INT)
    m_ij = L.MultiIndex([i, j], [p, q])

    body = [L.AssignAdd(y[i], A[m_ij] * x[j])]
    body = [L.ForRange(i, 0, p, body=body)]
    code += [L.ForRange(j, 0, q, body=body)]

    # Format into C and compile with CFFI
    format = Formatter(dtype=dtype)
    c_scalar = dtype_to_c_type(dtype)
    decl = f"void gemm({c_scalar} *y, {c_scalar} *A, {c_scalar} *x)"
    c_code = decl + "{\n" + format(L.StatementList(code)) + "\n}\n"

    ffibuilder = FFI()
    ffibuilder.cdef(decl + ";")
    ffibuilder.set_source(f"_gemv_{c_scalar}", c_code)
    ffibuilder.compile(verbose=True)
    _gemv = importlib.import_module(f"_gemv_{c_scalar}")
    gemv = _gemv.lib.gemm
    ffi = _gemv.ffi

    y = np.arange(p, dtype=dtype)
    x = np.arange(q, dtype=dtype)
    A = np.outer(y, x)

    py = ffi.cast(f"{c_scalar} *", y.ctypes.data)
    pA = ffi.cast(f"{c_scalar} *", A.ctypes.data)
    px = ffi.cast(f"{c_scalar} *", x.ctypes.data)

    # Compute expected result
    s2 = q * (q - 1) * (2 * q - 1) // 6 + 1
    result = np.arange(p, dtype=dtype) * s2

    gemv(py, pA, px)
    assert np.all(y == result)
