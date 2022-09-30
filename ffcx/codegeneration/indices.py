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
    lang: types.ModuleType
    ranges: typing.List[int]
    dim: int
    strides: typing.List[int]
    indices: typing.List

    def __init__(self, language, indices, ranges):
        self.lang = language
        self.ranges = ranges
        self.indices = indices
        self.dim = len(self.ranges)

        # Improve computation of stride
        self.strides = numpy.ones(self.dim, dtype=int)
        for i in range(self.dim - 1):
            self.strides[i] = numpy.prod(self.ranges[i + 1:])

        if 0 in self.ranges:
            self.strides[:] = 0

    def global_idx(self):
        if self.dim == 1:
            return self.indices[0]
        else:
            global_factors = [self.strides[i] * self.indices[i] for i in range(self.dim) if self.ranges[i] != 1]
            return self.lang.Sum(global_factors)

    def local_idx(self, idx):
        assert idx < self.dim
        return self.indices[idx]

    def global_size(self):
        return numpy.prod(self.ranges)

    def intersection(self, other):
        indices = []
        ranges = []
        for (idx, size) in zip(self.indices, self.ranges):
            if idx in other.indices:
                indices.append(idx)
                ranges.append(size)
        return MultiIndex(self.lang, indices, ranges)

    def union(self, other):
        indices = self.indices.copy()
        ranges = self.ranges.copy()
        for (idx, size) in zip(other.indices, other.ranges):
            if idx not in indices:
                indices.append(idx)
                ranges.append(size)
        return MultiIndex(self.lang, indices, ranges)

    def difference(self, other):
        indices = []
        ranges = []
        for (idx, size) in zip(self.indices, self.ranges):
            if idx not in other.indices:
                indices.append(idx)
                ranges.append(size)
        return MultiIndex(self.lang, indices, ranges)

    def collapse(self, name):
        indices = [self.lang.Symbol(name)]
        ranges = [numpy.prod(self.ranges, dtype=int)]
        return MultiIndex(self.lang, indices, ranges)

    def __hash__(self):
        return hash(self.global_idx())


def create_quadrature_index(lang: types.ModuleType, quadrature_rule: QuadratureRule) -> MultiIndex:
    ranges = [0]
    name = "iq"
    indices = [lang.Symbol(name)]
    if quadrature_rule:
        ranges = [quadrature_rule.weights.size]
        if quadrature_rule.has_tensor_factors:
            dim = len(quadrature_rule.tensor_factors)
            ranges = [factor[1].size for factor in quadrature_rule.tensor_factors]
            indices = [lang.Symbol(name + f"{i}") for i in range(dim)]

    iq = MultiIndex(lang, indices, ranges)
    return iq


def create_dof_index(lang: types.ModuleType, table_ref: UniqueTableReference, name: str) -> MultiIndex:
    if table_ref.has_tensor_factorisation:
        dim = len(table_ref.tensor_factors)
        ranges = [factor.values.shape[-1] for factor in table_ref.tensor_factors]
        indices = [lang.Symbol(name + f"{i}") for i in range(dim)]
    else:
        ranges = [table_ref.values.shape[-1]]
        indices = [lang.Symbol(name)]

    ic = MultiIndex(lang, indices, ranges)
    return ic
