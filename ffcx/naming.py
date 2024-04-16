# Copyright (C) 2009-2020 Anders Logg and Michal Habera
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Naming."""

from __future__ import annotations

import hashlib

import basix.ufl
import numpy as np
import numpy.typing as npt
import ufl

import ffcx
import ffcx.codegeneration


def compute_signature(
    ufl_objects: list[
        ufl.Form | basix.ufl._ElementBase | tuple[ufl.core.expr.Expr, npt.NDArray[np.float64]]
    ],
    tag: str,
) -> str:
    """Compute the signature hash.

    Based on the UFL type of the objects and an additional optional 'tag'.
    """
    object_signature = ""
    for ufl_object in ufl_objects:
        # Get signature from ufl object
        if isinstance(ufl_object, ufl.Form):
            kind = "form"
            object_signature += ufl_object.signature()
        elif isinstance(ufl_object, ufl.AbstractFiniteElement):
            object_signature += repr(ufl_object)
            kind = "element"
        elif isinstance(ufl_object, tuple) and isinstance(ufl_object[0], ufl.core.expr.Expr):
            expr = ufl_object[0]
            points = ufl_object[1]

            # FIXME Move this to UFL, cache the computation
            coeffs = ufl.algorithms.extract_coefficients(expr)
            consts = ufl.algorithms.analysis.extract_constants(expr)
            args = ufl.algorithms.analysis.extract_arguments(expr)

            rn = dict()
            rn.update(dict((c, i) for i, c in enumerate(coeffs)))
            rn.update(dict((c, i) for i, c in enumerate(consts)))
            rn.update(dict((c, i) for i, c in enumerate(args)))

            domains: list[ufl.Mesh] = []
            for coeff in coeffs:
                domains.append(*coeff.ufl_domains())
            for arg in args:
                domains.append(*arg.ufl_function_space().ufl_domains())
            for gc in ufl.algorithms.analysis.extract_type(expr, ufl.classes.GeometricQuantity):
                domains.append(*gc.ufl_domains())
            for const in consts:
                domains.append(const.ufl_domain())
            domains = ufl.algorithms.analysis.unique_tuple(domains)
            rn.update(dict((d, i) for i, d in enumerate(domains)))

            # Hash on UFL signature and points
            signature = ufl.algorithms.signature.compute_expression_signature(expr, rn)
            object_signature += signature
            object_signature += repr(points)

            kind = "expression"
        else:
            raise RuntimeError(f"Unknown ufl object type {ufl_object.__class__.__name__}")

    # Build combined signature
    signatures = [
        object_signature,
        str(ffcx.__version__),
        ffcx.codegeneration.get_signature(),
        kind,
        tag,
    ]
    string = ";".join(signatures)
    return hashlib.sha1(string.encode("utf-8")).hexdigest()


def integral_name(
    original_form: ufl.form.Form,
    integral_type: str,
    form_id: int,
    subdomain_id: tuple[int, ...] | tuple[str],
    prefix: str,
) -> str:
    """Get integral name."""
    sig = compute_signature([original_form], str((prefix, integral_type, form_id, subdomain_id)))
    return f"integral_{sig}"


def form_name(original_form: ufl.form.Form, form_id: int, prefix: str) -> str:
    """Get form name."""
    sig = compute_signature([original_form], str((prefix, form_id)))
    return f"form_{sig}"


def finite_element_name(ufl_element: basix.ufl._ElementBase, prefix: str) -> str:
    """Get finite element name."""
    assert isinstance(ufl_element, basix.ufl._ElementBase)
    sig = compute_signature([ufl_element], prefix)
    return f"element_{sig}"


def dofmap_name(ufl_element: basix.ufl._ElementBase, prefix: str) -> str:
    """Get DOF map name."""
    assert isinstance(ufl_element, basix.ufl._ElementBase)
    sig = compute_signature([ufl_element], prefix)
    return f"dofmap_{sig}"


def expression_name(
    expression: tuple[ufl.core.expr.Expr, npt.NDArray[np.floating]], prefix: str
) -> str:
    """Get expression name."""
    assert isinstance(expression[0], ufl.core.expr.Expr)
    sig = compute_signature([expression], prefix)
    return f"expression_{sig}"
