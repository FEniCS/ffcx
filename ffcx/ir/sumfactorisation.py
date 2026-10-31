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
import itertools
import math
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


def is_point_independent(factor: UniqueTableReferenceT) -> bool:
    """Check whether a tensor factor table has the same values at every point."""
    return bool(np.all(factor.values == factor.values[:, :, :1, :]))


def evaluation_order(
    factors: list[UniqueTableReferenceT], num_points: list[int]
) -> tuple[list[int], list[int]]:
    """Return how a linear combination of dofs is evaluated on a tensor product rule.

    Returns the directions in the order their dofs are contracted, and the
    number of points stored in each direction. A derivative of a low-order
    basis is often constant in its own direction, so that direction is
    stored at one point. Point-independent directions are contracted
    first: forming their differences of dofs before interpolating the
    others leaves fewer values to interpolate. The others follow, last
    direction first.
    """
    points = [1 if is_point_independent(f) else n for f, n in zip(factors, num_points)]
    directions = sorted(reversed(range(len(points))), key=lambda d: points[d] != 1)
    return directions, points


def evaluation_stage_sizes(
    num_dofs: list[int], num_points: list[int], directions: list[int]
) -> list[int]:
    """Return the sizes of the stages evaluating a linear combination of dofs, in order.

    Each stage holds, in each direction, the dofs or, once that direction is
    contracted, its stored points.
    """
    shape, sizes = list(num_dofs), []
    for d in directions:
        shape[d] = num_points[d]
        sizes.append(math.prod(shape))
    return sizes


def contraction_stage_size(direction: int, num_points: list[int]) -> int:
    """Return the size of a stage contracting an integrand factor over directions <= direction.

    The stage holds the points of the later directions.
    """
    return math.prod(num_points[direction + 1 :])


def _names(tables: list[UniqueTableReferenceT], direction: int) -> tuple[str, ...]:
    """Return the names of the 1D tables in a direction."""
    return tuple(_factors(td)[direction].name for td in tables)


def _stage_keys(factor: int, tables: list[UniqueTableReferenceT]) -> list[tuple]:
    """Identify an integrand factor contracted over each number of directions.

    Key d identifies the factor contracted with the tables over directions
    < d. The last key, over all directions, identifies the contribution to
    the element tensor.
    """
    keys: list[tuple] = [(factor,)]
    for direction in range(len(_factors(tables[0]))):
        keys.append(keys[-1] + (_names(tables, direction),))
    return keys


class ContractionTerm(typing.NamedTuple):
    """An element tensor block and the integrand factor contracted into it."""

    factor: int  # Index of the integrand factor in the factorisation graph
    blockmap: tuple[tuple[int, ...], ...]
    blockdata: BlockDataT
    tables: list[UniqueTableReferenceT]  # The argument tables
    stage_keys: list[tuple]  # See `_stage_keys`


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

    # Batches of groups. A batch caches its integrand factors on the whole
    # rule and contracts them, so batching bounds the space the caches take.
    batches: list[list[ContractionGroup]]
    # Products of the 1D quadrature weights and the 1D tables of the first
    # arguments, [points][dofs of each argument], by name
    tables: dict[str, npt.NDArray[np.float64]]
    # Names of the products, by direction and 1D table names
    table_names: dict[tuple[int, tuple[str, ...]], str]
    # 1D argument tables the contraction reads without weights
    factors: set[str]

    def weighted_table(self, direction: int, tables: list[UniqueTableReferenceT]) -> str:
        """Return the name of the weighted product of the 1D tables in a direction."""
        return self.table_names[(direction, _names(tables, direction))]


def _plan_group(terms: list[ContractionTerm]) -> ContractionGroup:
    """Plan the rows of blocks sharing their dof loops."""
    rank = len(terms[0].tables)

    # Contributions, and the contributions each block receives
    contributions: dict[tuple, ContractionTerm] = {}
    targets: dict[tuple, tuple[ContractionTerm, list[tuple]]] = {}
    for term in terms:
        key = term.stage_keys[-1]
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
            else:
                single += members
        rows.append(ContractionRow(row_targets, combined, single))
    return ContractionGroup(terms, rows)


def _group_size(group: ContractionGroup, num_points: list[int]) -> int:
    """Return the number of scalars the arrays of a group's contraction take.

    See `IntegralGenerator._sum_factorised_contraction`: a stage array per
    distinct stage, and per row an `h` array per combined set of
    contributions and an `acc` array.
    """
    last = len(num_points) - 1
    stages = {key for term in group.terms for key in term.stage_keys[1 : last + 1]}
    size = sum(contraction_stage_size(len(key) - 2, num_points) for key in stages)
    num_inner = int(_factors(group.terms[0].tables[-1])[last].values.shape[-1])
    return size + sum(
        len(row.combined) * num_points[last] * num_inner + num_inner for row in group.rows
    )


def plan_sum_factorisation(
    quadrature_rule: QuadratureRule,
    block_contributions: dict[tuple[tuple[int, ...], ...], list[BlockDataT]],
    max_scalars: int | None,
) -> SumFactorisationPlan | None:
    """Plan the contraction of the blocks of a sum-factorisable integrand.

    The integrand factors are cached and contracted in batches whose
    arrays take at most `max_scalars` scalars. Returns `None` if a single
    factor needs more.
    """
    dim = len(quadrature_rule.tensor_factors)
    last = dim - 1
    num_points = quadrature_rule.num_points

    # Group blocks with the same dof loops, so they share one loop nest
    grouped: dict[tuple, list[ContractionTerm]] = collections.defaultdict(list)
    for blockmap, blocks in sorted(block_contributions.items()):
        for blockdata in blocks:
            factor = blockdata.factor_indices_comp_indices[0][0]
            tables = [mad.tabledata for mad in blockdata.ma_data]
            term = ContractionTerm(factor, blockmap, blockdata, tables, _stage_keys(factor, tables))
            sizes = tuple(
                tuple(factor.values.shape[-1] for factor in _factors(td)) for td in term.tables
            )
            grouped[sizes].append(term)

    def plan_batch(batch: list[int]) -> list[ContractionGroup]:
        """Plan the groups contracting some of the integrand factors."""
        groups = []
        for terms in grouped.values():
            batch_terms = [term for term in terms if term.factor in batch]
            if batch_terms:
                groups.append(_plan_group(batch_terms))
        return groups

    def fits(batch: list[int]) -> bool:
        """Whether the caches of a batch of factors and the arrays contracting them fit."""
        if max_scalars is None:
            return True
        size = len(batch) * math.prod(num_points)
        return size + max(_group_size(g, num_points) for g in plan_batch(batch)) <= max_scalars

    # Add the factors, in order of use, to a batch while they fit
    batches: list[list[int]] = []
    for factor in dict.fromkeys(t.factor for terms in grouped.values() for t in terms):
        if batches and fits(batches[-1] + [factor]):
            batches[-1].append(factor)
        elif fits([factor]):
            batches.append([factor])
        else:
            return None
    planned = [plan_batch(batch) for batch in batches]

    # Directions with the same 1D rule share weighted tables
    table_names: dict[tuple[int, tuple[str, ...]], str] = {}
    by_values: dict[tuple[tuple[str, ...], bytes], str] = {}
    values: dict[str, npt.NDArray[np.float64]] = {}
    factors: set[str] = set()

    def weight(direction: int, tables: list[UniqueTableReferenceT]):
        """Add the product of 1D tables and the 1D quadrature weights in a direction."""
        names = _names(tables, direction)
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
        table_names[(direction, names)] = by_values[key]

    # The tables the contraction reads
    for group in itertools.chain.from_iterable(planned):
        for term in group.terms:
            for direction in range(last):
                weight(direction, term.tables)
        for row in group.rows:
            for members in row.combined:
                weight(last, members[0].tables[:1])
                factors.update(_factors(td)[last].name for td in members[0].tables[1:-1])
                factors.update(_factors(m.tables[-1])[last].name for m in members)
            for term in row.single:
                weight(last, term.tables)

    return SumFactorisationPlan(planned, values, table_names, factors)
