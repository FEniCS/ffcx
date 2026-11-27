# Copyright (C) 2026 Jørgen S. Dokken and Paul T. Kühner
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""External operators, evaluated in the kernels by the C functions they provide.

A function that maps the operands at the points pointwise gives the same
operator as interpolating the same map of the operands, which is the reference.
"""

import basix.ufl
import numpy as np
import pytest
import ufl

import ffcx.codegeneration.jit
import ffcx.compiler
import ffcx.naming
import ffcx.options
from ffcx.codegeneration.C.external import CExternalOperator

# UFL has no rule estimating the degree of an external operator yet, see the
# docstring of the handler in `ffcx.analysis`.
pytestmark = pytest.mark.filterwarnings("ignore:Missing degree estimation handler:UserWarning")


class CFunction(ufl.ExternalOperator):
    """An external operator evaluated by a C function, see `CExternalOperator`."""

    def __init__(
        self, *operands, function_space, name, code=None, derivatives=None, argument_slots=()
    ):
        """Initialise, with the name and the code of the C function."""
        self.c_name, self.c_code = name, code
        super().__init__(
            *operands,
            function_space=function_space,
            derivatives=derivatives,
            argument_slots=argument_slots,
        )

    def _ufl_expr_reconstruct_(self, *operands, add_kwargs={}, **kwargs):
        """Reconstruct, keeping the function."""
        add_kwargs = {"name": self.c_name, "code": self.c_code, **add_kwargs}
        return super()._ufl_expr_reconstruct_(*operands, add_kwargs=add_kwargs, **kwargs)

    # Defining __eq__ would otherwise drop the hash UFL computes and caches
    __hash__ = ufl.ExternalOperator.__hash__

    def _ufl_compute_hash_(self):
        """Hash, telling apart operators evaluated by different functions."""
        return hash((super()._ufl_compute_hash_(), self.c_name, self.c_code))

    def __eq__(self, other):
        """Check for equality, including the function."""
        return bool(
            super().__eq__(other)
            and isinstance(other, CFunction)
            and (self.c_name, self.c_code) == (other.c_name, other.c_code)
        )


SIGNATURE = """
void {name}(double* restrict values, const double* restrict operands, const double* restrict c,
            const double* restrict coordinate_dofs, const int* restrict entity_local_index,
            const uint8_t* restrict quadrature_permutation, void* custom_data)
{{
  // The interpolation points of the element of the operator
  for (int p = 0; p < {num_points}; ++p)
  {{
    {body}
  }}
}}
"""
# Into P2 and P3 on a triangle, 6 and 10 points
SQUARE = SIGNATURE.format(
    name="square", num_points=6, body="values[p] = operands[p] * operands[p];"
)
CUBE = SIGNATURE.format(
    name="cube", num_points=10, body="values[p] = operands[p] * operands[p] * operands[p];"
)
# The operands f and grad(g), three values at each point, to f grad(g), into vector P1
SCALE = SIGNATURE.format(
    name="scale",
    num_points=3,
    body="for (int i = 0; i < 2; ++i)\n"
    "      values[2 * p + i] = operands[3 * p] * operands[3 * p + 1 + i];",
)
CELLS = np.array(
    [[0.1, 0.2, 0], [2.0, 0.3, 0], [0.4, 1.7, 0], [2.1, 1.9, 0], [2.0, 0.3, 0], [0.4, 1.7, 0]]
)


def spaces():
    """A triangle mesh, a P2 space and the P1 and vector P1 spaces of the operators."""
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
    V = ufl.FunctionSpace(domain, basix.ufl.element("Lagrange", "triangle", 2))
    W = ufl.FunctionSpace(domain, basix.ufl.element("Lagrange", "triangle", 1))
    W2 = ufl.FunctionSpace(domain, basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
    return V, W, W2


def tabulate(module, form, integral_type, w, size):
    """Assemble the one integral of a form, on the first cell or the facet between both."""
    ffi = module.ffi
    integral = form.form_integrals[form.form_integral_offsets[getattr(module.lib, integral_type)]]
    A = np.zeros(size)
    coords = CELLS.flatten() if integral_type == "interior_facet" else CELLS[:3].flatten()
    entity, perm = np.zeros(2, dtype=np.int32), np.zeros(2, dtype=np.uint8)
    integral.tabulate_tensor_float64(
        ffi.cast("double *", A.ctypes.data),
        ffi.cast("double *", w.ctypes.data),
        ffi.NULL,
        ffi.cast("double *", coords.ctypes.data),
        ffi.cast("int *", entity.ctypes.data),
        ffi.cast("uint8_t *", perm.ctypes.data),
        ffi.NULL,
    )
    return A


def test_external_operator(compile_args):
    """External operators against the interpolations of the same maps, and the exact integral.

    Two functions in different spaces, and on both sides of an interior facet,
    each side built from its own cell. A vector-valued operator of two operands,
    ``(f, grad(g))``, with ``g`` numbered first, so that the function reads the
    operands in the order given, not in that of the coefficients. For a linear
    f, f**2 in P2 and f**3 in P3 are exact, which gives the integral of their
    sum exactly.
    """
    V, _, W2 = spaces()
    P3 = ufl.FunctionSpace(V.ufl_domain(), basix.ufl.element("Lagrange", "triangle", 3))
    g, f, v = ufl.Coefficient(V), ufl.Coefficient(V), ufl.TestFunction(V)
    square = CFunction(f, function_space=V, name="square", code=SQUARE)
    cube = CFunction(f, function_space=P3, name="cube", code=CUBE)
    scale = CFunction(f, ufl.grad(g), function_space=W2, name="scale", code=SCALE)
    assert isinstance(square, CExternalOperator)
    interpolate = ufl.Interpolate
    # UFL can't estimate the degree of an external operator: that of f**3 v
    dx, dS = ufl.dx(degree=5), ufl.dS(degree=5)
    pairs = [
        ((square + cube) * v * dx, (interpolate(f**2, V) + interpolate(f**3, P3)) * v * dx),
        (
            ufl.inner(scale, ufl.grad(v)) * dx,
            ufl.inner(interpolate(f * ufl.grad(g), W2), ufl.grad(v)) * dx,
        ),
        (
            (square("+") * v("+") + cube("-") * v("-")) * dS,
            (interpolate(f**2, V)("+") * v("+") + interpolate(f**3, P3)("-") * v("-")) * dS,
        ),
    ]
    forms = [form for pair in pairs for form in pair] + [(square + cube) * ufl.dx(degree=3)]
    compiled, module, _ = ffcx.codegeneration.jit.compile_forms(
        forms, cffi_extra_compile_args=compile_args
    )
    # Different on the two cells, so that reading the wrong side shows
    w = np.random.default_rng(0).uniform(0.5, 1.5, 12)
    for i, integral_type in enumerate(["cell", "cell", "interior_facet"]):
        size = 12 if integral_type == "interior_facet" else 6
        A = tabulate(module, compiled[2 * i], integral_type, w, size)
        A_ref = tabulate(module, compiled[2 * i + 1], integral_type, w, size)
        assert np.abs(A_ref).max() > 1e-3
        np.testing.assert_allclose(A, A_ref, rtol=1e-12)

    # The dofs of a linear f, its values at the points of P2 on the first cell
    x0, J = CELLS[0, :2], (CELLS[1:3, :2] - CELLS[0, :2]).T
    gradient = np.array([0.5, -0.7])
    w = 0.3 + (x0 + V.ufl_element().basix_element.points @ J.T) @ gradient
    # Exact for the cubic f**3 on the affine cell
    points, weights = basix.make_quadrature(basix.CellType.triangle, 3)
    values = 0.3 + (x0 + points @ J.T) @ gradient
    exact = abs(np.linalg.det(J)) * weights @ (values**2 + values**3)
    np.testing.assert_allclose(tabulate(module, compiled[-1], "cell", w, 1), [exact], rtol=1e-12)


def test_external_operator_code():
    """A function without code is declared, and the function and space are in the signature.

    The external operator is not a coefficient of the form. UFL leaves the
    function and the space out of its signature, so FFCx adds them.
    """
    V, W, _ = spaces()
    f, v = ufl.Coefficient(V), ufl.TestFunction(V)
    form = CFunction(f, function_space=W, name="square") * v * ufl.dx
    code, _ = ffcx.compiler.compile_ufl_objects([form], options=ffcx.options.get_options())
    assert "void square(" in "".join(code)
    analysis = ffcx.analysis.analyze_ufl_objects([form], np.float64)
    assert analysis.form_data[0].original_form.coefficients() == (f,)

    def signature(name, space, code=None):
        form = CFunction(f, function_space=space, name=name, code=code) * v * ufl.dx
        return ffcx.naming.compute_signature([form], "")

    signatures = {
        signature("square", W),
        signature("cube", W),
        signature("square", V),
        signature("square", W, SQUARE),
        signature("square", W, CUBE),
    }
    assert len(signatures) == 5


def test_external_operator_unsupported():
    """A Piola space, an operator without a C function, and one function defined twice."""
    V, W, _ = spaces()
    f, v = ufl.Coefficient(V), ufl.TestFunction(V)
    N1curl = ufl.FunctionSpace(W.ufl_domain(), basix.ufl.element("N1curl", "triangle", 1))
    cases = [
        (
            [
                ufl.inner(CFunction(ufl.grad(f), function_space=N1curl, name="g"), ufl.grad(v))
                * ufl.dx
            ],
            NotImplementedError,
            "identity pull back",
        ),
        (
            [ufl.ExternalOperator(f, function_space=W) * v * ufl.dx],
            NotImplementedError,
            "does not provide",
        ),
        (
            [
                CFunction(f, function_space=W, name="square", code=c) * v * ufl.dx
                for c in (SQUARE, CUBE)
            ],
            ValueError,
            "different code",
        ),
    ]
    for forms, error, match in cases:
        with pytest.raises(error, match=match):
            ffcx.compiler.compile_ufl_objects(forms, options=ffcx.options.get_options())
