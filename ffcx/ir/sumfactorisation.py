# Copyright (C) 2026 Garth N. Wells and Jørgen S. Dokken
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Plan the sum-factorised contraction of a cell integral on a tensor product rule.

The integrand factors are cached on the whole tensor product rule, then
contracted with the 1D argument tables one direction at a time. The plan
decides which element tensor blocks share loops and sums, and so which 1D
tables the contraction reads.
"""

from __future__ import annotations

import collections
import typing

import numpy as np
import numpy.typing as npt

from ffcx.ir.elementtables import UniqueTableReferenceT
from ffcx.ir.representationutils import QuadratureRule

if typing.TYPE_CHECKING:
    from ffcx.ir.integral import BlockDataT


def factorises(td: UniqueTableReferenceT, dim: int) -> bool:
    """Return whether a table factorises over the directions of a tensor product rule."""
    return (
        td.ttype not in ("zeros", "ones")
        and td.tensor_factors is not None
        and len(td.tensor_factors) == dim
    )


def _factors(td: UniqueTableReferenceT) -> list[UniqueTableReferenceT]:
    """Return the 1D tables of a table that factorises."""
    assert td.tensor_factors is not None
    return td.tensor_factors


def tables_key(tables: list[UniqueTableReferenceT], direction: int) -> tuple:
    """Identify the 1D tables of all arguments in a direction in a contraction.

    The first argument's tables have the 1D quadrature weights folded in.
    """
    return tuple((_factors(td)[direction].name, arg == 0) for arg, td in enumerate(tables))


class ContractionTerm(typing.NamedTuple):
    """An element tensor block and the integrand factor contracted into it."""

    factor: int  # Index of the integrand factor in the factorisation graph
    blockmap: tuple[tuple[int, ...], ...]
    blockdata: BlockDataT

    @property
    def tables(self) -> list[UniqueTableReferenceT]:
        """The argument tables."""
        return [mad.tabledata for mad in self.blockdata.ma_data]

    def stage_keys(self) -> list[tuple]:
        """Identify the integrand factor contracted over each number of directions.

        Key d identifies the factor contracted with the term's tables over
        directions < d. The last key, over all directions, identifies the
        term's contribution to the element tensor.
        """
        keys: list[tuple] = [(self.factor,)]
        for direction in range(len(_factors(self.tables[0]))):
            keys.append(keys[-1] + (tables_key(self.tables, direction),))
        return keys


class ContractionRow(typing.NamedTuple):
    """Element tensor blocks receiving the same contributions, which share one sum."""

    targets: list[ContractionTerm]  # A term of each block
    # Contributions sharing the 1D tables of all but the last argument in
    # the last direction, first combined over the last argument's table
    combined: list[list[ContractionTerm]]
    single: list[ContractionTerm]


class ContractionGroup(typing.NamedTuple):
    """Element tensor blocks sharing their dof loops."""

    terms: list[ContractionTerm]
    rows: list[ContractionRow]


class SumFactorisationPlan(typing.NamedTuple):
    """The contraction of a cell integral on a tensor product rule."""

    groups: list[ContractionGroup]
    # Products of the 1D quadrature weights and the 1D tables of the first
    # arguments, by direction and table names
    weighted_tables: dict[tuple[int, tuple[str, ...]], str]
    # 1D argument tables the contraction reads without weights
    factors: set[str]

    def weighted_table(self, direction: int, tables: list[UniqueTableReferenceT]) -> str:
        """Return the name of the weighted product of the 1D tables in a direction."""
        names = tuple(_factors(td)[direction].name for td in tables)
        return self.weighted_tables[(direction, names)]


def plan_sum_factorisation(
    quadrature_rule: QuadratureRule,
    block_contributions: dict[tuple[tuple[int, ...], ...], list[BlockDataT]],
) -> tuple[SumFactorisationPlan, dict[str, npt.NDArray[np.float64]]]:
    """Plan the contraction of the blocks of a sum-factorisable integrand.

    Returns the plan and the weighted tables it reads.
    """
    dim = len(quadrature_rule.tensor_factors)
    last = dim - 1

    # Group blocks with the same dof loops, so they share one loop nest
    grouped: dict[tuple, list[ContractionTerm]] = collections.defaultdict(list)
    for blockmap, blocks in sorted(block_contributions.items()):
        for blockdata in blocks:
            term = ContractionTerm(blockdata.factor_indices_comp_indices[0][0], blockmap, blockdata)
            sizes = tuple(
                tuple(factor.values.shape[-1] for factor in _factors(td)) for td in term.tables
            )
            grouped[sizes].append(term)

    # Directions with the same 1D rule share weighted tables
    weighted_tables: dict[tuple[int, tuple[str, ...]], str] = {}
    by_values: dict[tuple[tuple[str, ...], bytes], str] = {}
    values: dict[str, npt.NDArray[np.float64]] = {}
    factors: set[str] = set()

    def weight(direction: int, tables: list[UniqueTableReferenceT]):
        """Add the product of 1D tables and the 1D quadrature weights in a direction."""
        names = tuple(_factors(td)[direction].name for td in tables)
        weights = quadrature_rule.tensor_factors[direction][1]
        key = (names, weights.tobytes())
        if key not in by_values:
            table = weights.reshape((-1,) + (1,) * len(tables))
            for arg, td in enumerate(tables):
                factor = _factors(td)[direction].values
                shape = [weights.size] + [1] * len(tables)
                shape[arg + 1] = factor.shape[-1]
                table = table * factor.reshape(shape)
            name = f"FE_TP{len(values)}_Q{quadrature_rule.id()}"
            values[name] = table
            by_values[key] = name
        weighted_tables[(direction, names)] = by_values[key]

    groups = []
    for terms in grouped.values():
        rank = len(terms[0].tables)
        for term in terms:
            for direction in range(last):
                weight(direction, term.tables)

        # Contributions, and the contributions each block receives
        contributions: dict[tuple, ContractionTerm] = {}
        targets: dict[tuple, tuple[ContractionTerm, list[tuple]]] = {}
        for term in terms:
            key = term.stage_keys()[-1]
            contributions.setdefault(key, term)
            targets.setdefault(term.blockmap, (term, []))[1].append(key)

        signatures: dict[tuple, tuple[list[tuple], list[ContractionTerm]]] = {}
        for target, keys in targets.values():
            signatures.setdefault(tuple(sorted(keys, key=repr)), (keys, []))[1].append(target)

        rows = []
        for keys, row_targets in signatures.values():
            by_outer_tables: dict[tuple, list[ContractionTerm]] = collections.defaultdict(list)
            for key in keys:
                by_outer_tables[key[-1][:-1]].append(contributions[key])
            combined, single = [], []
            for members in by_outer_tables.values():
                if rank > 1 and len(members) > 1:
                    combined.append(members)
                    weight(last, members[0].tables[:1])
                    factors.update(_factors(td)[last].name for td in members[0].tables[1:-1])
                    factors.update(_factors(m.tables[-1])[last].name for m in members)
                else:
                    single += members
                    for member in members:
                        weight(last, member.tables)
            rows.append(ContractionRow(row_targets, combined, single))
        groups.append(ContractionGroup(terms, rows))

    return SumFactorisationPlan(groups, weighted_tables, factors), values
