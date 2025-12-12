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

from ffcx.ir.representation import FormIR


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
    names = []
    ids = []
    offsets = [0]
    domains = []
    # Note: the order of this list is defined by the enum ufcx_integral_type in ufcx.h
    for itg_type in ("cell", "exterior_facet", "interior_facet", "vertex", "ridge"):
        unsorted_integrals = []
        unsorted_ids = []
        unsorted_domains = []
        for name, _domains, id in zip(
            ir.integral_names[itg_type],
            ir.integral_domains[itg_type],
            ir.subdomain_ids[itg_type],
        ):
            unsorted_integrals += [name]
            unsorted_ids += [id]
            unsorted_domains += [_domains]

        id_sort = np.argsort(unsorted_ids)
        names += [unsorted_integrals[i] for i in id_sort]
        ids += [unsorted_ids[i] for i in id_sort]
        domains += [unsorted_domains[i] for i in id_sort]

        offsets.append(sum(len(d) for d in domains))

    return IntegralData(names, ids, offsets, domains)
