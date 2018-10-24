# -*- coding: utf-8 -*-
# Copyright (C) 2011-2017 Martin Sandve Alnæs
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Algorithms for factorizing argument dependent monomials."""

import logging
from functools import singledispatch

from ffc import FFCError
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
        fi = len(F.e2i)
        F.add_node(fi, expression=expr)
        F.e2i[expr] = fi
    return fi


# Reuse these empty objects where appropriate to save memory
noargs = {}


@singledispatch
def handler(v, fac, sf, F):
    # Error checking
    if any(fac):
        raise FFCError(
            "Assuming that a {0} cannot be applied to arguments. If this is wrong please report a bug.".
            format(type(v)))
    # Record non-argument subexpression
    raise FFCError("No arguments")


@handler.register(Sum)
def handle_sum(v, fac, sf, F):
    if len(fac) != 2:
        raise FFCError("Assuming binary sum here. This can be fixed if needed.")

    fac0 = fac[0]
    fac1 = fac[1]
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
                f0 = F.nodes[fi0]['expression']
                f1 = F.nodes[fi1]['expression']
                fisum = add_to_fv(f0 + f1, F)
            factors[argkey] = fisum

    else:  # non-arg + non-arg
        raise FFCError("No arguments")

    return factors


@handler.register(Product)
def handle_product(v, fac, sf, F):
    if len(fac) != 2:
        raise FFCError("Assuming binary product here. This can be fixed if needed.")
    fac0 = fac[0]
    fac1 = fac[1]

    if not fac0 and not fac1:  # non-arg * non-arg
        raise FFCError("No arguments")

    elif not fac0:  # non-arg * arg
        # Record products of non-arg operand with each factor of arg-dependent operand
        f0 = sf[0]
        factors = {}
        for k1 in sorted(fac1):
            fi1 = fac1[k1]
            factors[k1] = add_to_fv(f0 * F.nodes[fi1]['expression'], F)

    elif not fac1:  # arg * non-arg
        # Record products of non-arg operand with each factor of arg-dependent operand
        f1 = sf[1]
        factors = {}
        for k0 in sorted(fac0):
            f0 = F.nodes[fac0[k0]]['expression']
            factors[k0] = add_to_fv(f1 * f0, F)

    else:  # arg * arg
        # Record products of each factor of arg-dependent operand
        factors = {}
        for k0 in sorted(fac0):
            f0 = F.nodes[fac0[k0]]['expression']
            for k1 in sorted(fac1):
                f1 = F.nodes[fac1[k1]]['expression']
                argkey = tuple(sorted(k0 + k1))  # sort key for canonical representation
                factors[argkey] = add_to_fv(f0 * f1, F)

    return factors


@handler.register(Conj)
def handle_conj(v, fac, sf, F):

    fac = fac[0]
    if fac:
        factors = {}
        for k in fac:
            f0 = F.nodes[fac[k]]['expression']
            factors[k] = add_to_fv(Conj(f0), F)
    else:
        raise FFCError("No arguments")

    return factors


@handler.register(Division)
def handle_division(v, fac, sf, F):
    fac0 = fac[0]
    fac1 = fac[1]
    assert not fac1, "Cannot divide by arguments."

    if fac0:  # arg / non-arg
        # Record products of non-arg operand with each factor of arg-dependent operand
        f1 = sf[1]
        factors = {}
        for k0 in sorted(fac0):
            f0 = F.nodes[fac0[k0]]['expression']
            factors[k0] = add_to_fv(f0 / f1, F)

    else:  # non-arg / non-arg
        raise FFCError("No arguments")

    return factors


@handler.register(Conditional)
def handle_conditional(v, fac, sf, F):
    fac0 = fac[0]
    fac1 = fac[1]
    fac2 = fac[2]
    assert not fac0, "Cannot have argument in condition."

    if not (fac1 or fac2):  # non-arg ? non-arg : non-arg
        raise FFCError("No arguments")
    else:
        f0 = sf[0]
        f1 = sf[1]
        f2 = sf[2]

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
            f1 = z if fi1 is None else F.nodes[fi1]['expression']
            f2 = z if fi2 is None else F.nodes[fi2]['expression']
            factors[k] = add_to_fv(conditional(f0, f1, f2), F)

    return factors


def compute_argument_factorization(S, rank):
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

    # Factorize each subexpression in order:
    for si, attr in S.nodes.items():
        deps = S.out_edges[si]
        v = attr['expression']

        if si in arg_indices:
            assert len(deps) == 0
            # v is a modified Argument
            factors = {(si, ): one_index}
        else:
            fac = [S.nodes[d]['factors'] for d in deps]
            if not any(fac):
                # Entirely scalar (i.e. no arg factors)
                # Just add unchanged to F
                add_to_fv(v, F)
                factors = noargs
            else:
                # Get scalar factors for dependencies
                # which do not have arg factors
                sf = []
                for i, d in enumerate(deps):
                    if fac[i]:
                        sf.append(None)
                    else:
                        sf.append(S.nodes[d]['expression'])
                # Use appropriate handler to deal with Sum, Product, etc.
                factors = handler(v, fac, sf, F)

        attr['factors'] = factors

    assert len(F.nodes) == len(F.e2i)

    # Find the (only) node in S that is marked as 'target'
    # Should always be the last one.
    S_targets = [i for i, v in S.nodes.items() if v.get('target', False)]
    assert len(S_targets) == 1
    S_target = S_targets[0]

    # Get the factorizations of the target values
    if S.nodes[S_target]['factors'] == {}:
        if rank == 0:
            # Functionals and expressions: store as no args * factor
            factors = {(): F.e2i[S.nodes[S_target]['expression']]}
        else:
            # Zero form of arity 1 or higher: make factors empty
            factors = {}
    else:
        # Forms of arity 1 or higher:
        # Map argkeys from indices into SV to indices into AV,
        # and resort keys for canonical representation
        factors = {
            tuple(sorted(arg_indices.index(si) for si in argkey)): fi
            for argkey, fi in S.nodes[S_target]['factors'].items()
        }
    # Expecting all term keys to have length == rank
    # (this assumption will eventually have to change if we
    # implement joint bilinear+linear form factorization here)
    assert all(len(k) == rank for k in factors)

    # Indices into F that are needed for final result
    for i in factors.values():
        F.nodes[i]['target'] = []
    for k in factors:
        i = factors[k]
        F.nodes[i]['target'] += [k]

    # Compute dependencies in FV
    for i, v in F.nodes.items():
        expr = v['expression']
        if not expr._ufl_is_terminal_ and not expr._ufl_is_terminal_modifier_:
            for o in expr.ufl_operands:
                F.add_edge(i, F.e2i[o])

    return AV, F
