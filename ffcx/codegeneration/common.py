# Copyright (C) 2025 Paul T. Kühner
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

"""Common code for backend implemenations."""

import string
from typing import NamedTuple

import basix
import numpy as np

from ffcx.ir.representation import FormIR, IntegralIR


def template_keys(template: str) -> set[str]:
    """Set of expected data keys of a template."""
    return set(fname for _, fname, _, _ in string.Formatter().parse(template) if fname)


class IntegralData(NamedTuple):
    """Sorted integral data."""

    names: list[str]
    ids: list[int]
    offsets: list[int]
    domains: list[list[basix.CellType]]


def integral_data(ir: FormIR) -> IntegralData:
    """Extracts sorted intergral data from form."""
    names, ids, domains = [], [], []
    offsets = [0]
    # Note: the order of this list is defined by the enum ufcx_integral_type in ufcx.h
    for itg_type in ("cell", "exterior_facet", "interior_facet", "vertex", "ridge"):
        _ids = ir.subdomain_ids[itg_type]
        id_sort = np.argsort(_ids)

        ids += [_ids[i] for i in id_sort]
        names += [ir.integral_names[itg_type][i] for i in id_sort]
        domains += [ir.integral_domains[itg_type][i] for i in id_sort]

        offsets.append(offsets[-1] + sum(len(d) for d in domains[offsets[-1] :]))

    return IntegralData(names, ids, offsets, domains)


class KernelTensorSizes(NamedTuple):
    """Size information of kernel input data."""
    A: int
    w: int
    c: int
    coords: int
    local_index: int
    permutation: int

def integral_tensor_sizes(ir: IntegralIR) -> KernelTensorSizes:
    """Compute tensor sizes of integral IR input data."""
    A = np.prod(ir.expression.tensor_shape)
    w = sum(coeff.ufl_element().dim for coeff in ir.expression.coefficient_offsets.keys())
    c = sum(
        np.prod(constant.ufl_shape, dtype=int)
        for constant in ir.expression.original_constant_offsets.keys()
    )
    coords = ir.expression.number_coordinate_dofs * 3
    local_index = 2  # TODO: this is just an upper bound
    permutation = 2 if ir.expression.needs_facet_permutations else 0

    return KernelTensorSizes(A, w, c, coords, local_index, permutation)
