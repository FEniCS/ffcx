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
# along with UFLACS. If not, see <http://www.gnu.org/licenses/>

"""Rebuilding UFL expressions from linearized representation of computational graph."""

import numpy

from ufl import product
from ufl.permutation import compute_indices

from ufl import as_vector
from ufl.classes import MultiIndex, IndexSum, Product
from ufl.corealg.multifunction import MultiFunction
from ufl.utils.indexflattening import flatten_multiindex, shape_to_strides

from ffc.log import error
from ffc.uflacs.analysis.modified_terminals import is_modified_terminal


class ReconstructScalarSubexpressions(MultiFunction):

    def __init__(self):
        super(ReconstructScalarSubexpressions, self).__init__()

    # No fallbacks, need to specify each type or group of types explicitly
    def expr(self, o, *args, **kwargs):
        error("No handler for type %s" % type(o))

    def terminal(self, o):
        error("Not expecting terminal expression in here, got %s." % type(o))

    # These types are not expected to be part of the graph at this point
    def unexpected(self, o, *args, **kwargs):
        error("Not expecting expression of type %s in here." % type(o))
    multi_index = unexpected
    expr_list = unexpected
    expr_mapping = unexpected
    utility_type = unexpected
    label = unexpected
    component_tensor = unexpected
    list_tensor = unexpected
    transposed = unexpected
    variable = unexpected

    def scalar_nary(self, o, ops):
        if o.ufl_shape != ():
            error("Expecting scalar.")
        sops = [op[0] for op in ops]
        return [o._ufl_expr_reconstruct_(*sops)]

    # Unary scalar functions
    math_function = scalar_nary
    abs = scalar_nary
    min_value = scalar_nary
    max_value = scalar_nary
    # Binary scalar functions
    power = scalar_nary
    bessel_function = scalar_nary  # TODO: Is this ok?
    atan_2 = scalar_nary

    def condition(self, o, ops):
        # A condition is always scalar, so len(op) == 1
        sops = [op[0] for op in ops]
        return [o._ufl_expr_reconstruct_(*sops)]

    def conditional(self, o, ops):
        # A condition can be non scalar
        symbols = []
        n = len(ops[1])
        if len(ops[0]) != 1:
            error("Condition should be scalar.")
        if n != len(ops[2]):
            error("Conditional branches should have same shape.")
        for i in range(len(ops[1])):
            sops = (ops[0][0], ops[1][i], ops[2][i])
            symbols.append(o._ufl_expr_reconstruct_(*sops))
        return symbols

    def division(self, o, ops):
        if len(ops) != 2:
            error("Expecting two operands.")
        if len(ops[1]) != 1:
            error("Expecting scalar divisor.")
        b, = ops[1]
        return [o._ufl_expr_reconstruct_(a, b) for a in ops[0]]

    def sum(self, o, ops):
        if len(ops) != 2:
            error("Expecting two operands.")
        if len(ops[0]) != len(ops[1]):
            error("Expecting scalar divisor.")
        return [o._ufl_expr_reconstruct_(a, b) for a, b in zip(ops[0], ops[1])]

    def product(self, o, ops):
        if len(ops) != 2:
            error("Expecting two operands.")

        # Get the simple cases out of the way
        if len(ops[0]) == 1:  # True scalar * something
            a, = ops[0]
            return [Product(a, b) for b in ops[1]]
        elif len(ops[1]) == 1:  # Something * true scalar
            b, = ops[1]
            return [Product(a, b) for a in ops[0]]

        # Neither of operands are true scalars, this is the tricky part
        o0, o1 = o.ufl_operands

        # Get shapes and index shapes
        fi = o.ufl_free_indices
        fi0 = o0.ufl_free_indices
        fi1 = o1.ufl_free_indices
        fid = o.ufl_index_dimensions
        fid0 = o0.ufl_index_dimensions
        fid1 = o1.ufl_index_dimensions

        # Need to map each return component to one component of o0 and one component of o1.
        indices = compute_indices(fid)

        # Compute which component of o0 is used in component (comp,ind) of o
        # Compute strides within free index spaces
        ist0 = shape_to_strides(fid0)
        ist1 = shape_to_strides(fid1)
        # Map o0 and o1 indices to o indices
        indmap0 = [fi.index(i) for i in fi0]
        indmap1 = [fi.index(i) for i in fi1]
        indks = [(flatten_multiindex([ind[i] for i in indmap0], ist0),
                  flatten_multiindex([ind[i] for i in indmap1], ist1))
                 for ind in indices]

        # Build products for scalar components
        results = [Product(ops[0][k0], ops[1][k1]) for k0, k1 in indks]
        return results

    def index_sum(self, o, ops):
        summand, mi = o.ufl_operands
        ic = mi[0].count()
        fi = summand.ufl_free_indices
        fid = summand.ufl_index_dimensions
        ipos = fi.index(ic)
        d = fid[ipos]

        # Compute "macro-dimensions" before and after i in the total shape of a
        predim = product(summand.ufl_shape) * product(fid[:ipos])
        postdim = product(fid[ipos+1:])

        # Map each flattened total component of summand to
        # flattened total component of indexsum o by removing
        # axis corresponding to summation index ii.
        ss = ops[0]  # Scalar subexpressions of summand
        if len(ss) != predim * postdim * d:
            error("Mismatching number of subexpressions.")
        sops = []
        for i in range(predim):
            iind = i * (postdim * d)
            for k in range(postdim):
                ind = iind + k
                sops.append([ss[ind + j * postdim] for j in range(d)])

        # For each scalar output component, sum over collected subcomponents
        # TODO: Need to split this into binary additions to work with future CRSArray format,
        #       i.e. emitting more expressions than there are symbols for this node.
        results = [sum(sop) for sop in sops]
        return results


    # TODO: To implement compound tensor operators such as dot and inner,
    # we need to identify which index to do the contractions over,
    # and build expressions such as sum(a*b for a,b in zip(aops, bops))


def rebuild_expression_from_graph(G):
    "This is currently only used by tests."
    w = rebuild_with_scalar_subexpressions(G)

    # Find expressions of final v
    if len(w) == 1:
        return w[0]
    else:
        return as_vector(w)  # TODO: Consider shape of initial v


def rebuild_with_scalar_subexpressions(G, targets=None):
    """Build a new expression2index mapping where each subexpression is scalar valued.

    Input:
    - G.e2i
    - G.V
    - G.V_symbols
    - G.total_unique_symbols

    Output:
    - NV   - Array with reverse mapping from index to expression
    - nvs  - Tuple of ne2i indices corresponding to the last vertex of G.V

    Old output now no longer returned but possible to restore if needed:
    - ne2i - Mapping from scalar subexpressions to a contiguous unique index
    - W    - Array with reconstructed scalar subexpressions for each original symbol
    """

    # From simplefsi3d.ufl:
    # GRAPH SIZE: len(G.V), G.total_unique_symbols
    # GRAPH SIZE: 16251   635272
    # GRAPH SIZE:   473     8210
    # GRAPH SIZE:  9663   238021
    # GRAPH SIZE: 88913  3448634  #  3.5 M!!!

    # Algorithm to apply to each subexpression
    reconstruct_scalar_subexpressions = ReconstructScalarSubexpressions()

    # Array to store the scalar subexpression in for each symbol
    W = numpy.empty(G.total_unique_symbols, dtype=object)

    # Iterate over each graph node in order
    for i, v in enumerate(G.V):
        # Find symbols of v components
        vs = G.V_symbols[i]

        # Skip if there's nothing new here (should be the case for indexing types)
        if all(W[s] is not None for s in vs):
            continue

        if is_modified_terminal(v):
            # if v.ufl_free_indices:
            #     error("Expecting no free indices.")
            sh = v.ufl_shape
            if sh:
                # Store each terminal expression component.
                # We may not actually need all of these later,
                # but that will be optimized away.
                # Note: symmetries will be dealt with in the value numbering.
                ws = [v[c] for c in compute_indices(sh)]
            else:
                # Store single modified terminal expression component
                if len(vs) != 1:
                    error("Expecting single symbol for scalar valued modified terminal.")
                ws = [v]
            # FIXME: Replace ws[:] with 0's if its table is empty
            # Possible redesign: loop over modified terminals only first,
            # then build tables for them, set W[s] = 0.0 for modified terminals with zero table,
            # then loop over non-(modified terminal)s to reconstruct expression.
        else:
            # Find symbols of operands
            sops = []
            for j, vop in enumerate(v.ufl_operands):
                if isinstance(vop, MultiIndex):
                    # TODO: Store MultiIndex in G.V and allocate a symbol to it for this to work
                    if not isinstance(v, IndexSum):
                        error("Not expecting a %s." % type(v))
                    sops.append(())
                else:
                    # TODO: Build edge datastructure and use instead?
                    # k = G.E[i][j]
                    k = G.e2i[vop]
                    sops.append(G.V_symbols[k])

            # Fetch reconstructed operand expressions
            wops = [tuple(W[k] for k in so) for so in sops]

            # Reconstruct scalar subexpressions of v
            ws = reconstruct_scalar_subexpressions(v, wops)

            # Store all scalar subexpressions for v symbols
            if len(vs) != len(ws):
                error("Expecting one symbol for each expression.")

        # Store each new scalar subexpression in W at the index of its symbol
        handled = set()
        for s, w in zip(vs, ws):
            if W[s] is None:
                W[s] = w
                handled.add(s)
            else:
                assert s in handled  # Result of symmetry!

    # Find symbols of requested targets or final v from input graph
    if targets is None:
        targets = [G.V[-1]]

    # Attempt to extend this to multiple target expressions
    scalar_target_expressions = []
    for target in targets:
        ti = G.e2i[target]
        vs = G.V_symbols[ti]
        # Sanity check: assert that we've handled these symbols
        if any(W[s] is None for s in vs):
            error("Expecting that all symbols in vs are handled at this point.")
        scalar_target_expressions.append([W[s] for s in vs])

    # Return the scalar expressions for each of the components
    assert len(scalar_target_expressions) == 1  # TODO: Currently expected by callers, fix those first
    return scalar_target_expressions[0]  # ... TODO: then return list
