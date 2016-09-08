# -*- coding: utf-8 -*-
# Copyright (C) 2011-2015 Martin Sandve Alnæs
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
# along with UFLACS. If not, see <http://www.gnu.org/licenses/>

"""Algorithms for factorizing argument dependent monomials."""

import numpy

from six import itervalues, iterkeys, iteritems
from six.moves import xrange as range

from ufl import as_ufl, conditional
from ufl.classes import Argument
from ufl.classes import Division
from ufl.classes import Product
from ufl.classes import Sum
from ufl.classes import Conditional
from ufl.classes import Zero
from ufl.algorithms import extract_type

from ffc.log import ffc_assert, error

from uflacs.analysis.graph_dependencies import compute_dependencies
from uflacs.analysis.modified_terminals import analyse_modified_terminal, strip_modified_terminal


def _build_arg_sets(V):
    "Build arg_sets = { argument number: set(j for j where V[j] is a modified Argument with this number) }"
    arg_sets = {}
    for i, v in enumerate(V):
        arg = strip_modified_terminal(v)
        if not isinstance(arg, Argument):
            continue
        num = arg.number()
        arg_set = arg_sets.get(num)
        if arg_set is None:
            arg_set = {}
            arg_sets[num] = arg_set
        arg_set[i] = v
    return arg_sets


def _build_argument_indices_from_arg_sets(V, arg_sets):
    "Build ordered list of indices to modified arguments."
    # Build set of all indices of V referring to modified arguments
    arg_indices = set()
    for js in itervalues(arg_sets):
        arg_indices.update(js)

    # Make a canonical ordering of vertex indices for modified arguments
    def arg_ordering_key(i):
        "Return a key for sorting argument vertex indices based on the properties of the modified terminal."
        mt = analyse_modified_terminal(arg_ordering_key.V[i])
        arg = mt.terminal
        assert isinstance(arg, Argument)
        assert arg.number() >= 0
        return (arg.number(),
                arg.part(),
                mt.reference_value,
                mt.component,
                mt.global_derivatives,
                mt.local_derivatives,
                mt.restriction,
                mt.averaged)
    arg_ordering_key.V = V
    ordered_arg_indices = sorted(arg_indices, key=arg_ordering_key)

    return ordered_arg_indices


def build_argument_indices(V):
    "Build ordered list of indices to modified arguments."
    arg_sets = _build_arg_sets(V)
    ordered_arg_indices = _build_argument_indices_from_arg_sets(V, arg_sets)
    return ordered_arg_indices


def build_argument_dependencies(dependencies, arg_indices):
    "Preliminary algorithm: build list of argument vertex indices each vertex (indirectly) depends on."
    n = len(dependencies)
    A = [[] for i in range(n)]  # TODO: Use array
    for i, deps in enumerate(dependencies):
        argdeps = []
        for j in deps:
            if j in arg_indices:
                argdeps.append(j)
            else:
                argdeps.extend(A[j])
        A[i] = sorted(argdeps)
    return A


class Factors(object): # TODO: Refactor code in this file by using a class like this
    def __init__(self):
        self.FV = []
        self.e2fi = {}

    def add(self, expr):
        add_to_fv(expr, self.FV, self.e2fi)


def add_to_fv(expr, FV, e2fi):
    "Add expression expr to factor vector FV and expr->FVindex mapping e2fi."
    fi = e2fi.get(expr)
    if fi is None:
        fi = len(e2fi)
        FV.append(expr)
        e2fi[expr] = fi
    return fi


# Reuse these empty objects where appropriate to save memory
noargs = {}


def handle_modified_terminal(i, v, F, FV, e2fi, arg_indices, AV, sv2av):
    # v is a modified terminal...
    if i in arg_indices:
        # ... a modified Argument
        argkey = (i,)
        fi = None

        # Adding 1 as an expression allows avoiding special representation by representing "v" as "1*v"
        one = add_to_fv(as_ufl(1.0), FV, e2fi)
        factors = {argkey: one}

        assert AV[sv2av[i]] == v
    else:
        # ... record a non-argument modified terminal
        factors = noargs
        fi = add_to_fv(v, FV, e2fi)
    return fi, factors


def handle_sum(i, v, deps, F, FV, sv2fv, e2fi):
    ffc_assert(len(deps) == 2, "Assuming binary sum here. This can be fixed if needed.")
    fac0 = F[deps[0]]
    fac1 = F[deps[1]]

    argkeys = sorted(set(iterkeys(fac0)) | set(iterkeys(fac1)))

    if argkeys:  # f*arg + g*arg = (f+g)*arg
        keylen = len(argkeys[0])
        fi = None
        factors = {}
        for argkey in argkeys:
            ffc_assert(len(argkey) == keylen, "Expecting equal argument rank terms among summands.")

            fi0 = fac0.get(argkey)
            fi1 = fac1.get(argkey)
            if fi0 is None:
                fisum = fi1
            elif fi1 is None:
                fisum = fi0
            else:
                f0 = FV[fi0]
                f1 = FV[fi1]
                fisum = add_to_fv(f0 + f1, FV, e2fi)
            factors[argkey] = fisum

    else:  # non-arg + non-arg
        factors = noargs
        fi = add_to_fv(v, FV, e2fi)

    return fi, factors


def handle_product(i, v, deps, F, FV, sv2fv, e2fi):
    ffc_assert(len(deps) == 2, "Assuming binary product here. This can be fixed if needed.")
    fac0 = F[deps[0]]
    fac1 = F[deps[1]]

    if not fac0 and not fac1:  # non-arg * non-arg
        # Record non-argument product
        factors = noargs
        f0 = FV[sv2fv[deps[0]]]
        f1 = FV[sv2fv[deps[1]]]
        assert f1 * f0 == v
        fi = add_to_fv(v, FV, e2fi)
        assert FV[fi] == v

    elif not fac0:  # non-arg * arg
        # Record products of non-arg operand with each factor of arg-dependent operand
        f0 = FV[sv2fv[deps[0]]]
        factors = {}
        for k1 in sorted(fac1):
            fi1 = fac1[k1]
            factors[k1] = add_to_fv(f0 * FV[fi1], FV, e2fi)
        fi = None

    elif not fac1:  # arg * non-arg
        # Record products of non-arg operand with each factor of arg-dependent operand
        f1 = FV[sv2fv[deps[1]]]
        factors = {}
        for k0 in sorted(fac0):
            f0 = FV[fac0[k0]]
            factors[k0] = add_to_fv(f1 * f0, FV, e2fi)
        fi = None

    else:  # arg * arg
        # Record products of each factor of arg-dependent operand
        factors = {}
        for k0 in sorted(fac0):
            f0 = FV[fac0[k0]]
            for k1 in sorted(fac1):
                f1 = FV[fac1[k1]]
                argkey = tuple(sorted(k0 + k1))  # sort key for canonical representation
                factors[argkey] = add_to_fv(f0 * f1, FV, e2fi)
        fi = None

    return fi, factors


def handle_division(i, v, deps, F, FV, sv2fv, e2fi):
    fac0 = F[deps[0]]
    fac1 = F[deps[1]]
    assert not fac1, "Cannot divide by arguments."

    if fac0:  # arg / non-arg
        # Record products of non-arg operand with each factor of arg-dependent operand
        f1 = FV[sv2fv[deps[1]]]
        factors = {}
        for k0 in sorted(fac0):
            f0 = FV[fac0[k0]]
            factors[k0] = add_to_fv(f0 / f1, FV, e2fi)
        fi = None

    else:  # non-arg / non-arg
        # Record non-argument subexpression
        factors = noargs
        fi = add_to_fv(v, FV, e2fi)

    return fi, factors


def handle_conditional(i, v, deps, F, FV, sv2fv, e2fi):
    fac0 = F[deps[0]]
    fac1 = F[deps[1]]
    fac2 = F[deps[2]]
    assert not fac0, "Cannot have argument in condition."

    if not (fac1 or fac2):  # non-arg ? non-arg : non-arg
        # Record non-argument subexpression
        fi = add_to_fv(v, FV, e2fi)
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

        fi = None
        factors = {}

        z = as_ufl(0.0)
        zfi = add_to_fv(z, FV, e2fi)  # TODO: flake8 complains zfi is unused, is that ok?

        # In general, can decompose like this:
        #    conditional(c, sum_i fi*ui, sum_j fj*uj) -> sum_i conditional(c, fi, 0)*ui + sum_j conditional(c, 0, fj)*uj
        mas = sorted(set(fac1.keys()) | set(fac2.keys()))
        for k in mas:
            fi1 = fac1.get(k)
            fi2 = fac2.get(k)
            f1 = z if fi1 is None else FV[fi1]
            f2 = z if fi2 is None else FV[fi2]
            factors[k] = add_to_fv(conditional(f0, f1, f2), FV, e2fi)

    return fi, factors


def handle_operator(i, v, deps, F, FV, sv2fv, e2fi):
    # Error checking
    if any(F[deps[j]] for j in range(len(deps))):
        error("Assuming that a {0} cannot be applied to arguments. If this is wrong please report a bug.".format(type(v)))
    # Record non-argument subexpression
    fi = add_to_fv(v, FV, e2fi)
    factors = noargs
    return fi, factors


def collect_argument_factors(SV, dependencies, arg_indices):
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

          SV[-1] := sum(FV[fik] * product(AV[j] for j in aik) for aik, fik in IM.items())

        where := means equivalence in the mathematical sense,
        of course in a different technical representation.

    """
    # Extract argument component subgraph
    AV = [SV[j] for j in arg_indices]
    #av2sv = arg_indices
    sv2av = dict((j, i) for i, j in enumerate(arg_indices))
    assert all(AV[i] == SV[j] for i, j in enumerate(arg_indices))
    assert all(AV[i] == SV[j] for j, i in iteritems(sv2av))

    # Data structure for building non-argument factors
    FV = []
    e2fi = {}

    # Hack to later build dependencies for the FV entries that change K*K -> K**2
    two = add_to_fv(as_ufl(2), FV, e2fi)  # FIXME: Might need something more robust here

    # Intermediate factorization for each vertex in SV on the format
    # F[i] = None # if SV[i] does not depend on arguments
    # F[i] = { argkey: fi } # if SV[i] does depend on arguments, where:
    #   FV[fi] is the expression SV[i] with arguments factored out
    #   argkey is a tuple with indices into SV for each of the argument components SV[i] depends on
    # F[i] = { argkey1: fi1, argkey2: fi2, ... } # if SV[i] is a linear combination of multiple argkey configurations
    F = numpy.empty(len(SV), dtype=object)
    sv2fv = numpy.zeros(len(SV), dtype=int)

    # Factorize each subexpression in order:
    for i, v in enumerate(SV):
        deps = dependencies[i]

        # These handlers insert values in sv2fv and F
        if not len(deps):
            fi, factors = handle_modified_terminal(i, v, F, FV, e2fi, arg_indices, AV, sv2av)
        elif isinstance(v, Sum):
            fi, factors = handle_sum(i, v, deps, F, FV, sv2fv, e2fi)
        elif isinstance(v, Product):
            fi, factors = handle_product(i, v, deps, F, FV, sv2fv, e2fi)
        elif isinstance(v, Division):
            fi, factors = handle_division(i, v, deps, F, FV, sv2fv, e2fi)
        elif isinstance(v, Conditional):
            fi, factors = handle_conditional(i, v, deps, F, FV, sv2fv, e2fi)
        else:  # All other operators
            fi, factors = handle_operator(i, v, deps, F, FV, sv2fv, e2fi)

        if fi is not None:
            sv2fv[i] = fi
        F[i] = factors

    assert not noargs, "This dict was not supposed to be filled with anything!"

    # Throw away superfluous items in array
    # FV = FV[:len(e2fi)]
    assert len(FV) == len(e2fi)

    # Get the factorization of the final value # TODO: Support simultaneous factorization of multiple integrands?
    IM = F[-1]

    # Map argkeys from indices into SV to indices into AV, and resort keys for canonical representation
    IM = dict((tuple(sorted(sv2av[j] for j in argkey)), fi) for argkey, fi in iteritems(IM))

    # If this is a non-argument expression, point to the expression from IM (not sure if this is useful)
    if any([not AV, not IM, not arg_indices]):
        assert all([not AV, not IM, not arg_indices])
        IM = {(): len(FV) - 1}

    return FV, e2fi, AV, IM


def rebuild_scalar_graph_from_factorization(AV, FV, IM):
    # TODO: What about multiple target_variables?

    # Build initial graph
    SV = []
    SV.extend(AV)
    SV.extend(FV)
    se2i = dict((s, i) for i, s in enumerate(SV))

    def add_vertex(h):
        # Avoid adding vertices twice
        i = se2i.get(h)
        if i is None:
            se2i[h] = len(SV)
            SV.append(h)

    # Add factorization monomials
    argkeys = sorted(iterkeys(IM))
    fs = []
    for argkey in argkeys:
        # Start with coefficients
        f = FV[IM[argkey]]
        # f = 1

        # Add binary products with each argument in order
        for argindex in argkey:
            f = f * AV[argindex]
            add_vertex(f)

        # Add product with coefficients last
        # f = f*FV[IM[argkey]]
        # add_vertex(f)

        # f is now the full monomial, store it as a term for sum below
        fs.append(f)

    # Add sum of factorization monomials
    g = 0
    for f in fs:
        g = g + f
        add_vertex(g)

    # Rebuild dependencies
    dependencies = compute_dependencies(se2i, SV)

    return SV, se2i, dependencies


def compute_argument_factorization(SV, target_variables, dependencies):

    # TODO: Use target_variables! Currently just assuming the last vertex is the target here...

    if list(target_variables) != [len(SV) - 1]:
        ffc_assert(not extract_type(SV[-1], Argument),
                   "Multiple or nonscalar Argument dependent expressions not supported in factorization.")
        AV = []
        FV = SV
        IM = {}
        return AV, FV, IM, target_variables, dependencies

    assert list(target_variables) == [len(SV) - 1]

    arg_indices = build_argument_indices(SV)
    #A = build_argument_dependencies(dependencies, arg_indices)
    FV, e2fi, AV, IM = collect_argument_factors(SV, dependencies, arg_indices)

    # Indices into FV that are needed for final result
    target_variables = sorted(itervalues(IM))

    dependencies = compute_dependencies(e2fi, FV)

    return IM, AV, FV, target_variables, dependencies
