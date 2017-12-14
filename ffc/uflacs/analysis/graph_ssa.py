# -*- coding: utf-8 -*-
# Copyright (C) 2011-2017 Martin Sandve Alnæs
#
# This file is part of UFLACS.
#
# UFLACS is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# UFLACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with UFLACS. If not, see <http://www.gnu.org/licenses/>.

"""Algorithms for working with computational graphs."""

import numpy

from ufl.classes import (GeometricQuantity, ConstantValue,
                         Argument, Coefficient,
                         Grad, Restricted, Indexed,
                         MathFunction)
from ufl.checks import is_cellwise_constant
from ffc.log import error

from ffc.uflacs.analysis.crsarray import CRSArray


def default_partition_seed(expr, rank):
    """
    Partition 0: Piecewise constant on each cell (including Real and DG0 coefficients)
    Partition 1: Depends on x
    Partition 2: Depends on x and coefficients
    Partitions [3,3+rank): depend on argument with count partition-3
    """
    # TODO: Use named constants for the partition numbers here

    modifiers = (Grad, Restricted, Indexed)  # FIXME: Add CellAvg, FacetAvg types here, others?
    if isinstance(expr, modifiers):
        return default_partition_seed(expr.ufl_operands[0], rank)

    elif isinstance(expr, Argument):
        ac = expr.number()
        assert 0 <= ac < rank
        poffset = 3
        p = poffset + ac
        return p

    elif isinstance(expr, Coefficient):
        if is_cellwise_constant(expr):  # This is crap, doesn't include grad modifier
            return 0
        else:
            return 2

    elif isinstance(expr, GeometricQuantity):
        if is_cellwise_constant(expr):  # This is crap, doesn't include grad modifier
            return 0
        else:
            return 1

    elif isinstance(expr, ConstantValue):
        return 0

    else:
        error("Don't know how to handle %s" % expr)


def mark_partitions(V, active, dependencies, rank,
                    partition_seed=default_partition_seed,
                    partition_combiner=max):
    """FIXME: Cover this with tests.

    Input:
    - V            - Array of expressions.
    - active       - Boolish array.
    - dependencies - CRSArray with V dependencies.
    - partition_seed - Policy for determining the partition of a terminalish.
    - partition_combiner - Policy for determinging the partition of an operator.

    Output:
    - partitions   - Array of partition int ids.
    """
    n = len(V)
    assert len(active) == n
    assert len(dependencies) == n
    partitions = numpy.zeros(n, dtype=int)
    for i, v in enumerate(V):
        deps = dependencies[i]
        if active[i]:
            if len(deps):
                p = partition_combiner([partitions[d] for d in deps])
            else:
                p = partition_seed(v, rank)
        else:
            p = -1
        partitions[i] = p
    return partitions


"""
def build_factorized_partitions():
    num_points = [3]

    # dofrange = (begin, end)
    # dofblock = ()  |  (dofrange0,)  |  (dofrange0, dofrange1)

    partitions = {}

    # partitions["piecewise"] = partition of expressions independent of quadrature and argument loops
    partitions["piecewise"] = []

    # partitions["varying"][np] = partition of expressions dependent on np quadrature but independent of argument loops
    partitions["varying"] = dict((np, []) for np in num_points)

    # partitions["argument"][np][iarg][dofrange] = partition depending on this dofrange of argument iarg
    partitions["argument"] = dict((np, [dict() for i in range(rank)]) for np in num_points)

    # partitions["integrand"][np][dofrange] = partition depending on this dofrange of argument iarg
    partitions["integrand"] = dict((np, dict()) for np in num_points)
"""

def compute_dependency_count(dependencies):
    """FIXME: Test"""
    n = len(dependencies)
    depcount = numpy.zeros(n, dtype=int)
    for i in range(n):
        for d in dependencies[i]:
            depcount[d] += 1
    return depcount


def invert_dependencies(dependencies, depcount):
    """FIXME: Test"""
    n = len(dependencies)
    m = sum(depcount)
    invdeps = [()] * n
    for i in range(n):
        for d in dependencies[i]:
            invdeps[d] = invdeps[d] + (i,)
    return CRSArray.from_rows(invdeps, n, m, int)


def default_cache_score_policy(vtype, ndeps, ninvdeps, partition):
    # Start at 1 and then multiply with various heuristic factors
    s = 1

    # Is the type particularly expensive to compute?
    expensive = (MathFunction,)
    if vtype in expensive:  # Could make a type-to-cost mapping, but this should do.
        s *= 20

    # More deps roughly correlates to more operations
    s *= ndeps

    # If it is reused several times let that count significantly
    s *= ninvdeps ** 3  # 1->1, 2->8, 3->27

    # Finally let partition count for something?
    # Or perhaps we need some more information, such as
    # when x from outer loop is used by y within inner loop.

    return s


def compute_cache_scores(V, active, dependencies, inverse_dependencies, partitions,
                         cache_score_policy=default_cache_score_policy):
    """FIXME: Cover with tests.

    TODO: Experiment with heuristics later when we have functional code generation.
    """
    n = len(V)
    score = numpy.zeros(n, dtype=int)
    for i, v in enumerate(V):
        if active[i]:
            deps = dependencies[i]
            ndeps = len(deps)
            invdeps = inverse_dependencies[i]
            ninvdeps = len(invdeps)
            p = partitions[i]
            s = cache_score_policy(type(v), ndeps, ninvdeps, p)
        else:
            s = -1
        score[i] = s
    return score


import heapq


def allocate_registers(active, partitions, targets,
                       scores, max_registers, score_threshold):
    """FIXME: Cover with tests.

    TODO: Allow reuse of registers, reducing memory usage.

    TODO: Probably want to sort within partitions.
    """
    # Check for consistent number of variables
    n = len(scores)
    assert n == len(active)
    assert n == len(partitions)
    num_targets = len(targets)

    # Analyse scores
    #min_score = min(scores)
    #max_score = max(scores)
    #mean_score = sum(scores) // n

    # Can allocate a number of registers up to given threshold
    num_to_allocate = max(num_targets,
                          min(max_registers, n) - num_targets)
    to_allocate = set()

    # For now, just using an arbitrary heuristic algorithm to select m largest scores
    queue = [(-scores[i], i) for i in range(n) if active[i]]
    heapq.heapify(queue)

    # Always allocate registers for all targets, for simplicity
    # in the rest of the code generation pipeline
    to_allocate.update(targets)

    # Allocate one register each for max_registers largest symbols
    for r in range(num_to_allocate):
        s, i = heapq.heappop(queue)
        if -s <= score_threshold:
            break
        if i in targets:
            continue
        to_allocate.add(i)

    registers_used = len(to_allocate)

    # Some consistency checks
    assert num_to_allocate <= max_registers
    assert registers_used <= num_to_allocate + len(targets)
    assert registers_used <= max(max_registers, len(targets))

    # Mark allocations
    allocations = numpy.zeros(n, dtype=int)
    allocations[:] = -1
    for r, i in enumerate(sorted(to_allocate)):
        allocations[i] = r

    # Possible data structures for improved register allocations
    # register_status = numpy.zeros(max_registers, dtype=int)

    # Stack/set of free registers (should wrap in stack abstraction):
    # free_registers = numpy.zeros(max_registers, dtype=int)
    # num_free_registers = max_registers
    # free_registers[:] = reversed(xrange(max_registers))

    return allocations
