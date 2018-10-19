# -*- coding: utf-8 -*-
# Copyright (C) 2011-2017 Martin Sandve Alnæs
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Algorithms for factorizing argument dependent monomials."""

import logging
import itertools

import numpy

from ffc import FFCError
from ffc.uflacs.analysis.dependencies import compute_dependencies
from ffc.uflacs.analysis.graph import ExpressionGraph
from ffc.uflacs.analysis.modified_terminals import (analyse_modified_terminal,
                                                    strip_modified_terminal)
from ufl import as_ufl, conditional
from ufl.classes import Argument, Conditional, Division, Product, Sum, Zero, Conj

logger = logging.getLogger(__name__)


def build_argument_indices(S):
    """Build ordered list of indices to modified arguments."""

    arg_indices = []
    for i, v in S.nodes.items():
        arg = strip_modified_terminal(v['expression'])
        if isinstance(arg, Argument):
            arg_indices.append(i)

    # Make a canonical ordering of vertex indices for modified arguments
    def arg_ordering_key(i):
        """Return a key for sorting argument vertex indices based on
        the properties of the modified terminal."""
        mt = analyse_modified_terminal(S.nodes[i]['expression'])
        return mt.argument_ordering_key()

    ordered_arg_indices = sorted(arg_indices, key=arg_ordering_key)
    return ordered_arg_indices


def add_to_fv(expr, F):
    """Add expression expr to factor vector FV and expr->FVindex mapping e2fi."""
    fi = F.e2i.get(expr)
    if fi is None:
        fi = len(e2fi)
        F.add_node(fi, expression=expr)
        e2fi[expr] = fi
    return fi


# Reuse these empty objects where appropriate to save memory
noargs = {}


def handle_sum(v, si, deps, SV_factors, F, sv2fv):
    if len(deps) != 2:
        raise FFCError("Assuming binary sum here. This can be fixed if needed.")

    fac0 = SV_factors[deps[0]]
    fac1 = SV_factors[deps[1]]
    argkeys = set(fac0) | set(fac1)

    if argkeys:  # f*arg + g*arg = (f+g)*arg
        argkeys = sorted(argkeys)
        keylen = len(argkeys[0])
        factors = {}
        for argkey in argkeys:
            if len(argkey) != keylen:
                raise FFCError("Expecting equal argument rank terms among summands.")

            fi0 = fac0.get(argkey)
            fi1 = fac1.get(argkey)
            if fi0 is None:
                fisum = fi1
            elif fi1 is None:
                fisum = fi0
            else:
                f0 = FV[fi0]
                f1 = FV[fi1]
                fisum = add_to_fv(f0 + f1, F)
            factors[argkey] = fisum

    else:  # non-arg + non-arg
        factors = noargs
        sv2fv[si] = add_to_fv(v, F)

    return factors


def handle_product(v, si, deps, SV_factors, F, sv2fv):
    if len(deps) != 2:
        raise FFCError("Assuming binary product here. This can be fixed if needed.")
    fac0 = SV_factors[deps[0]]
    fac1 = SV_factors[deps[1]]

    if not fac0 and not fac1:  # non-arg * non-arg
        # Record non-argument product
        factors = noargs
        f0 = FV[sv2fv[deps[0]]]
        f1 = FV[sv2fv[deps[1]]]
        assert f1 * f0 == v
        sv2fv[si] = add_to_fv(v, F)
        assert FV[sv2fv[si]] == v

    elif not fac0:  # non-arg * arg
        # Record products of non-arg operand with each factor of arg-dependent operand
        f0 = FV[sv2fv[deps[0]]]
        factors = {}
        for k1 in sorted(fac1):
            fi1 = fac1[k1]
            factors[k1] = add_to_fv(f0 * FV[fi1], F)

    elif not fac1:  # arg * non-arg
        # Record products of non-arg operand with each factor of arg-dependent operand
        f1 = FV[sv2fv[deps[1]]]
        factors = {}
        for k0 in sorted(fac0):
            f0 = FV[fac0[k0]]
            factors[k0] = add_to_fv(f1 * f0, F)

    else:  # arg * arg
        # Record products of each factor of arg-dependent operand
        factors = {}
        for k0 in sorted(fac0):
            f0 = FV[fac0[k0]]
            for k1 in sorted(fac1):
                f1 = FV[fac1[k1]]
                argkey = tuple(sorted(k0 + k1))  # sort key for canonical representation
                factors[argkey] = add_to_fv(f0 * f1, F)

    return factors


def handle_conj(v, si, deps, SV_factors, F, sv2fv):

    fac = SV_factors[deps[0]]
    if fac:
        factors = {}
        for k in fac:
            f0 = FV[fac[k]]
            factors[k] = add_to_fv(Conj(f0), F)
    else:
        factors = noargs
        f0 = FV[sv2fv[deps[0]]]
        sv2fv[si] = add_to_fv(v, F)
        assert FV[sv2fv[si]] == v

    return factors


def handle_division(v, si, deps, SV_factors, F, sv2fv):
    fac0 = SV_factors[deps[0]]
    fac1 = SV_factors[deps[1]]
    assert not fac1, "Cannot divide by arguments."

    if fac0:  # arg / non-arg
        # Record products of non-arg operand with each factor of arg-dependent operand
        f1 = FV[sv2fv[deps[1]]]
        factors = {}
        for k0 in sorted(fac0):
            f0 = FV[fac0[k0]]
            factors[k0] = add_to_fv(f0 / f1, F)

    else:  # non-arg / non-arg
        # Record non-argument subexpression
        factors = noargs
        sv2fv[si] = add_to_fv(v, F)

    return factors


def handle_conditional(v, si, deps, SV_factors, F, sv2fv):
    fac0 = SV_factors[deps[0]]
    fac1 = SV_factors[deps[1]]
    fac2 = SV_factors[deps[2]]
    assert not fac0, "Cannot have argument in condition."

    if not (fac1 or fac2):  # non-arg ? non-arg : non-arg
        # Record non-argument subexpression
        sv2fv[si] = add_to_fv(v, F)
        factors = noargs
    else:
        f0 = FV[sv2fv[deps[0]]]
        f1 = FV[sv2fv[deps[1]]]
        f2 = FV[sv2fv[deps[2]]]

        # Term conditional(c, argument, non-argument) is not legal unless non-argument is 0.0
        assert fac1 or isinstance(f1, Zero)
        assert fac2 or isinstance(f2, Zero)
        assert () not in fac1
        assert () not in fac2

        z = as_ufl(0.0)

        # In general, can decompose like this:
        #    conditional(c, sum_i fi*ui, sum_j fj*uj) -> sum_i conditional(c, fi, 0)*ui + sum_j conditional(c, 0, fj)*uj
        mas = sorted(set(fac1.keys()) | set(fac2.keys()))
        factors = {}
        for k in mas:
            fi1 = fac1.get(k)
            fi2 = fac2.get(k)
            f1 = z if fi1 is None else FV[fi1]
            f2 = z if fi2 is None else FV[fi2]
            factors[k] = add_to_fv(conditional(f0, f1, f2), F)

    return factors


def handle_operator(v, si, deps, SV_factors, F, sv2fv):

    # Error checking
    if any(SV_factors[d] for d in deps):
        raise FFCError(
            "Assuming that a {0} cannot be applied to arguments. If this is wrong please report a bug.".
            format(type(v)))
    # Record non-argument subexpression
    sv2fv[si] = add_to_fv(v, F)
    factors = noargs
    return factors


def compute_argument_factorization(S, SV_target, rank):
    """Factorizes a scalar expression graph w.r.t. scalar Argument
    components.

    The result is a triplet (AV, FV, IM):

      - The scalar argument component subgraph:

          AV[ai] = v

        with the property

          SV[arg_indices] == AV[:]

      - An expression graph vertex list with all non-argument factors:

          FV[fi] = f

        with the property that none of the expressions depend on Arguments.

      - A dict representation of the final integrand of rank r:

          IM = { (ai1_1, ..., ai1_r): fi1, (ai2_1, ..., ai2_r): fi2, }

        This mapping represents the factorization of SV[-1] w.r.t. Arguments s.t.:

          SV[-1] := sum(FV[fik] * product(AV[ai] for ai in aik) for aik, fik in IM.items())

        where := means equivalence in the mathematical sense,
        of course in a different technical representation.

    """
    # Extract argument component subgraph
    arg_indices = build_argument_indices(S)
    AV = [S.nodes[i]['expression'] for i in arg_indices]
    sv2av = {si: ai for ai, si in enumerate(arg_indices)}

    # Data structure for building non-argument factors
    F = ExpressionGraph()

    # Adding 1.0 as an expression allows avoiding special representation
    # of arguments when first visited by representing "v" as "1*v"
    one_index = add_to_fv(as_ufl(1.0), F)

    # Intermediate factorization for each vertex in SV on the format
    # SV_factors[si] = None # if SV[si] does not depend on arguments
    # SV_factors[si] = { argkey: fi } # if SV[si] does depend on arguments, where:
    #   FV[fi] is the expression SV[si] with arguments factored out
    #   argkey is a tuple with indices into SV for each of the argument components SV[si] depends on
    # SV_factors[si] = { argkey1: fi1, argkey2: fi2, ... } # if SV[si]
    # is a linear combination of multiple argkey configurations
    SV_factors = numpy.empty(len(S.nodes), dtype=object)
    si2fi = numpy.zeros(len(S.nodes), dtype=int)

    # Factorize each subexpression in order:
    for si, attr in S.nodes.items():
        deps = S.out_edges[si]
        v = attr['expression']

        # These handlers insert values in si2fi and SV_factors
        if not len(deps):
            if si in arg_indices:
                # v is a modified Argument
                factors = {(si, ): one_index}
            else:
                # v is a modified non-Argument terminal
                si2fi[si] = add_to_fv(v, F)
                factors = noargs
        else:
            # These quantities could be better input args to handlers:
            # facs = [SV_factors[d] for d in deps]
            # fs = [FV[sv2fv[d]] for d in deps]
            if isinstance(v, Sum):
                handler = handle_sum
            elif isinstance(v, Conj):
                handler = handle_conj
            elif isinstance(v, Product):
                handler = handle_product
            elif isinstance(v, Division):
                handler = handle_division
            elif isinstance(v, Conditional):
                handler = handle_conditional
            else:  # All other operators
                handler = handle_operator
            factors = handler(v, si, deps, SV_factors, F, si2fi)

        SV_factors[si] = factors

    assert not noargs, "This dict was not supposed to be filled with anything!"

    assert len(FV) == len(e2fi)

    # Get the factorizations of the target values
    IMs = []
    if SV_factors[SV_target] == {}:
        if rank == 0:
            # Functionals and expressions: store as no args * factor
            factors = {(): si2fi[SV_target]}
        else:
            # Zero form of arity 1 or higher: make factors empty
            factors = {}
    else:
        # Forms of arity 1 or higher:
        # Map argkeys from indices into SV to indices into AV,
        # and resort keys for canonical representation
        factors = {
            tuple(sorted(sv2av[si] for si in argkey)): fi
            for argkey, fi in SV_factors[SV_target].items()
        }
    # Expecting all term keys to have length == rank
    # (this assumption will eventually have to change if we
    # implement joint bilinear+linear form factorization here)
    assert all(len(k) == rank for k in factors)
    IMs.append(factors)

    # Recompute dependencies in FV
    FV_deps = compute_dependencies(e2fi, FV)

    # Indices into FV that are needed for final result
    FV_targets = list(itertools.chain(sorted(IM.values()) for IM in IMs))

    return IMs, AV, FV, FV_deps, FV_targets
