# -*- coding: utf-8 -*-
# Copyright (C) 2009-2017 Anders Logg
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Some basics for generating C code."""

import hashlib

import ffc
import ufl


def make_name(prefix, basename, signature):
    pre = prefix.lower() + "_" if prefix else ""
    sig = str(signature).lower()
    return "{}{}_{}".format(pre, basename, sig)


def make_integral_name(prefix, integral_type, original_form, form_id, subdomain_id, parameters):
    sig = compute_signature([original_form], str(form_id), parameters)
    basename = "{}_integral_{}".format(integral_type, sig)
    return make_name(prefix, basename, subdomain_id)


def compute_signature(ufl_objects, tag, parameters, coordinate_mapping=False):
    """Compute the signature hash for jit modules.
    Based on the UFL type of the objects,
    relevant compilation parameters, signatures of this ffc version, ufc.h and ufc_geometry.h
    and an additional optional 'tag' (used with filename for command-line invocation).

    Note:
    ----
    The parameter `coordinate_mapping` is used to force compilation of finite element
    as a coordinate mapping element. There is no way to find this information
    just by looking at type of `ufl_object` passed.

    """

    object_signature = ""
    for ufl_object in ufl_objects:
        # Get signature from ufl object
        if isinstance(ufl_object, ufl.Form):
            kind = "form"
            object_signature += ufl_object.signature()
        elif isinstance(ufl_object, ufl.Mesh):
            # When coordinate mapping is represented by a Mesh, just getting
            # its coordinate element
            object_signature += repr(ufl_object.ufl_coordinate_element())
            kind = "coordinate_mapping"
        elif coordinate_mapping and isinstance(ufl_object, ufl.FiniteElementBase):
            object_signature += repr(ufl_object)
            kind = "coordinate_mapping"
        elif isinstance(ufl_object, ufl.FiniteElementBase):
            object_signature += repr(ufl_object)
            kind = "element"
        elif isinstance(ufl_object, tuple) and isinstance(ufl_object[0], ufl.core.expr.Expr):
            expr = ufl_object[0]

            # FIXME Move this to UFL, cache the computation
            coeffs = ufl.algorithms.extract_coefficients(expr)
            consts = ufl.algorithms.analysis.extract_constants(expr)

            rn = dict()
            rn.update(dict((c, i) for i, c in enumerate(coeffs)))
            rn.update(dict((c, i) for i, c in enumerate(consts)))

            domains = []
            for coeff in coeffs:
                domains.append(*coeff.ufl_domains())
            for gc in ufl.algorithms.analysis.extract_type(expr, ufl.classes.GeometricQuantity):
                domains.append(*gc.ufl_domains())

            domains = ufl.algorithms.analysis.unique_tuple(domains)
            rn.update(dict((d, i) for i, d in enumerate(domains)))

            siganture = ufl.algorithms.signature.compute_expression_signature(expr, rn)
            object_signature += siganture
            kind = "expression"
        else:
            raise RuntimeError("Unknown ufl object type {}".format(ufl_object.__class__.__name__))

    # Compute deterministic string of relevant parameters
    parameters_signature = ffc.parameters.compute_jit_signature(parameters)

    # Build combined signature
    signatures = [
        object_signature,
        parameters_signature,
        str(ffc.__version__),
        ffc.codegeneration.get_signature(),
        kind,
        tag
    ]
    string = ";".join(signatures)
    sig = hashlib.sha1(string.encode('utf-8')).hexdigest()

    return sig
