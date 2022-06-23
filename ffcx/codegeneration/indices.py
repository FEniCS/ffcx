# Copyright (C) 2022 Igor A. Baratta
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import types
import typing
import numpy

from ffcx.ir.elementtables import UniqueTableReference
from ffcx.ir.representationutils import QuadratureRule


class MultiIndex:
    L: types.ModuleType
    name: str
    ranges: typing.List[int]
    dim: int
    strides: typing.List[int]
    indices: typing.List

    def __init__(self, language, name, ranges):
        self.L = language
        self.name = name
        self.ranges = ranges
        self.dim = len(ranges)

        # Improve computation of stride
        self.strides = [numpy.prod(ranges[i + 1:], dtype=int) for i in range(self.dim - 1)]
        if 0 in self.ranges:
            self.strides.append(0)
        else:
            self.strides.append(1)

        if self.dim == 1:
            self.indices = [self.L.Symbol(name)]
        else:
            self.indices = [self.L.Symbol(name + f"{i}") for i in range(self.dim)]

    def global_idx(self):
        global_factors = [self.strides[i] * self.indices[i] for i in range(self.dim)]
        return self.L.Sum(global_factors)

    def local_idx(self, idx):
        assert idx < self.dim
        return self.indices[idx]


def create_quadrature_index(lang: types.ModuleType, quadrature_rule: QuadratureRule) -> MultiIndex:
    ranges = [0]
    if quadrature_rule:
        ranges = [quadrature_rule.weights.size]
        if quadrature_rule.has_tensor_factors:
            ranges = [factor[1].size for factor in quadrature_rule.tensor_factors]
    iq = MultiIndex(lang, "iq", ranges)
    return iq


def create_dof_index(lang: types.ModuleType, table_ref: UniqueTableReference, name: str) -> MultiIndex:
    ranges = [table_ref.values.shape[3]]
    if table_ref.has_tensor_factorisation:
        ranges = [factor.values.shape[-1] for factor in table_ref.tensor_factors]
    ic = MultiIndex(lang, name, ranges)
    return ic
