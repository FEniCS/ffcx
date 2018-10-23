# -*- coding: utf-8 -*-
# Copyright (C) 2011-2017 Martin Sandve Alnæs
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Linearized data structure for the computational graph."""

import logging

import numpy

import ufl
from ffc import FFCError
from ffc.uflacs.analysis.modified_terminals import is_modified_terminal
from ffc.uflacs.analysis.valuenumbering import ValueNumberer

logger = logging.getLogger(__name__)


class ExpressionGraph(object):
    """A directed multi-edge graph, allowing multiple edges
    between the same nodes, and respecting the insertion order
    of nodes and edges."""

    def __init__(self):

        # Data structures for directed multi-edge graph
        self.nodes = {}
        self.out_edges = {}
        self.in_edges = {}

        # Index to expression
        self.V = []

        # Expression to index dict
        self.e2i = {}

    def number_of_nodes(self):
        return len(self.nodes)

    def add_node(self, key, **kwargs):
        """Add a node with optional properties"""
        self.nodes[key] = kwargs
        self.out_edges[key] = []
        self.in_edges[key] = []

    def add_edge(self, node1, node2):
        """Add a directed edge from node1 to node2"""
        if node1 not in self.nodes or node2 not in self.nodes:
            raise KeyError("Adding edge to unknown node")

        self.out_edges[node1] += [node2]
        self.in_edges[node2] += [node1]


def build_graph_vertices(expression, scalar=False):
    # Count unique expression nodes

    G = ExpressionGraph()

    G.e2i = {}
    _count_nodes_with_unique_post_traversal(expression, G.e2i, scalar)

    # Invert the map to get index->expression
    G.V = sorted(G.e2i, key=G.e2i.get)

    # Add nodes to 'new' graph structure
    for i, v in enumerate(G.V):
        G.add_node(i, expression=v)

    # Get vertex index representing input expression root
    G.V_target = G.e2i[expression]

    return G


def build_scalar_graph(expression):
    """Build list representation of expression graph covering the given
    expressions.
    """

    # Populate with vertices
    G = build_graph_vertices(expression, scalar=False)

    # Build more fine grained computational graph of scalar subexpressions
    scalar_expression = rebuild_with_scalar_subexpressions(G)

    # Build new list representation of graph where all
    # vertices of V represent single scalar operations
    G = build_graph_vertices(scalar_expression, scalar=True)

    # Compute graph edges
    V_deps = []
    for v in G.V:
        if v._ufl_is_terminal_ or v._ufl_is_terminal_modifier_:
            V_deps.append(())
        else:
            V_deps.append([G.e2i[o] for o in v.ufl_operands])

    for i, edges in enumerate(V_deps):
        for j in edges:
            G.add_edge(i, j)

    return G


class ReconstructScalarSubexpressions(ufl.corealg.multifunction.MultiFunction):
    def __init__(self):
        super().__init__()

    # No fallbacks, need to specify each type or group of types explicitly
    def expr(self, o, *args, **kwargs):
        raise FFCError("No handler for type %s" % type(o))

    def terminal(self, o):
        raise FFCError("Not expecting terminal expression in here, got %s." % type(o))

    # These types are not expected to be part of the graph at this point
    def unexpected(self, o, *args, **kwargs):
        raise FFCError("Not expecting expression of type %s in here." % type(o))

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
            raise FFCError("Expecting scalar.")
        sops = [op[0] for op in ops]
        return [o._ufl_expr_reconstruct_(*sops)]

    # Unary scalar functions
    math_function = scalar_nary
    abs = scalar_nary
    min_value = scalar_nary
    max_value = scalar_nary
    real = scalar_nary
    imag = scalar_nary

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
            raise FFCError("Condition should be scalar.")
        if n != len(ops[2]):
            raise FFCError("Conditional branches should have same shape.")
        for i in range(len(ops[1])):
            sops = (ops[0][0], ops[1][i], ops[2][i])
            symbols.append(o._ufl_expr_reconstruct_(*sops))
        return symbols

    def conj(self, o, ops):
        if len(ops) != 1:
            raise FFCError("Expecting one operand")
        if o.ufl_shape != ():
            raise FFCError("Expecting scalar.")
        return [o._ufl_expr_reconstruct_(x) for x in ops[0]]

    def division(self, o, ops):
        if len(ops) != 2:
            raise FFCError("Expecting two operands.")
        if len(ops[1]) != 1:
            raise FFCError("Expecting scalar divisor.")
        b, = ops[1]
        return [o._ufl_expr_reconstruct_(a, b) for a in ops[0]]

    def sum(self, o, ops):
        if len(ops) != 2:
            raise FFCError("Expecting two operands.")
        if len(ops[0]) != len(ops[1]):
            raise FFCError("Expecting scalar divisor.")
        return [o._ufl_expr_reconstruct_(a, b) for a, b in zip(ops[0], ops[1])]

    def product(self, o, ops):
        if len(ops) != 2:
            raise FFCError("Expecting two operands.")

        # Get the simple cases out of the way
        if len(ops[0]) == 1:  # True scalar * something
            a, = ops[0]
            return [ufl.classes.Product(a, b) for b in ops[1]]
        elif len(ops[1]) == 1:  # Something * true scalar
            b, = ops[1]
            return [ufl.classes.Product(a, b) for a in ops[0]]

        # Neither of operands are true scalars, this is the tricky part
        o0, o1 = o.ufl_operands

        # Get shapes and index shapes
        fi = o.ufl_free_indices
        fi0 = o0.ufl_free_indices
        fi1 = o1.ufl_free_indices
        fid = o.ufl_index_dimensions
        fid0 = o0.ufl_index_dimensions
        fid1 = o1.ufl_index_dimensions

        # Need to map each return component to one component of o0 and
        # one component of o1
        indices = ufl.permutation.compute_indices(fid)

        # Compute which component of o0 is used in component (comp,ind) of o
        # Compute strides within free index spaces
        ist0 = ufl.utils.indexflattening.shape_to_strides(fid0)
        ist1 = ufl.utils.indexflattening.shape_to_strides(fid1)
        # Map o0 and o1 indices to o indices
        indmap0 = [fi.index(i) for i in fi0]
        indmap1 = [fi.index(i) for i in fi1]
        indks = [(ufl.utils.indexflattening.flatten_multiindex([ind[i] for i in indmap0], ist0),
                  ufl.utils.indexflattening.flatten_multiindex([ind[i] for i in indmap1], ist1))
                 for ind in indices]

        # Build products for scalar components
        results = [ufl.classes.Product(ops[0][k0], ops[1][k1]) for k0, k1 in indks]
        return results

    def index_sum(self, o, ops):
        summand, mi = o.ufl_operands
        ic = mi[0].count()
        fi = summand.ufl_free_indices
        fid = summand.ufl_index_dimensions
        ipos = fi.index(ic)
        d = fid[ipos]

        # Compute "macro-dimensions" before and after i in the total shape of a
        predim = ufl.product(summand.ufl_shape) * ufl.product(fid[:ipos])
        postdim = ufl.product(fid[ipos + 1:])

        # Map each flattened total component of summand to
        # flattened total component of indexsum o by removing
        # axis corresponding to summation index ii.
        ss = ops[0]  # Scalar subexpressions of summand
        if len(ss) != predim * postdim * d:
            raise FFCError("Mismatching number of subexpressions.")
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


def rebuild_with_scalar_subexpressions(G):
    """Build a new expression2index mapping where each subexpression is scalar valued.

    Input:
    - G.e2i
    - G.V
    - G.V_symbols
    - G.total_unique_symbols

    Output:
    - NV   - Array with reverse mapping from index to expression
    - nvs  - Tuple of ne2i indices corresponding to the last vertex of G.V
    """

    # Compute symbols over graph and rebuild scalar expression
    value_numberer = ValueNumberer(G)
    V_symbols = value_numberer.compute_symbols()
    total_unique_symbols = value_numberer.symbol_count

    # Algorithm to apply to each subexpression
    reconstruct_scalar_subexpressions = ReconstructScalarSubexpressions()

    # Array to store the scalar subexpression in for each symbol
    W = numpy.empty(total_unique_symbols, dtype=object)

    # Iterate over each graph node in order
    for i, v in enumerate(G.V):
        # Find symbols of v components
        vs = V_symbols[i]

        # Skip if there's nothing new here (should be the case for indexing types)
        if all(W[s] is not None for s in vs):
            continue

        if is_modified_terminal(v):
            # if v.ufl_free_indices:
            #     raise FFCError("Expecting no free indices.")
            sh = v.ufl_shape
            if sh:
                # Store each terminal expression component. We may not
                # actually need all of these later, but that will be
                # optimized away.
                # Note: symmetries will be dealt with in the value numbering.
                ws = [v[c] for c in ufl.permutation.compute_indices(sh)]
            else:
                # Store single modified terminal expression component
                if len(vs) != 1:
                    raise FFCError("Expecting single symbol for scalar valued modified terminal.")
                ws = [v]
            # FIXME: Replace ws[:] with 0's if its table is empty
            # Possible redesign: loop over modified terminals only first,
            # then build tables for them, set W[s] = 0.0 for modified terminals with zero table,
            # then loop over non-(modified terminal)s to reconstruct expression.
        else:
            # Find symbols of operands
            sops = []
            for j, vop in enumerate(v.ufl_operands):
                if isinstance(vop, ufl.classes.MultiIndex):
                    # TODO: Store MultiIndex in G.V and allocate a symbol to it for this to work
                    if not isinstance(v, ufl.classes.IndexSum):
                        raise FFCError("Not expecting a %s." % type(v))
                    sops.append(())
                else:
                    # TODO: Build edge datastructure and use instead?
                    # k = G.E[i][j]
                    k = G.e2i[vop]
                    sops.append(V_symbols[k])

            # Fetch reconstructed operand expressions
            wops = [tuple(W[k] for k in so) for so in sops]

            # Reconstruct scalar subexpressions of v
            ws = reconstruct_scalar_subexpressions(v, wops)

            # Store all scalar subexpressions for v symbols
            if len(vs) != len(ws):
                raise FFCError("Expecting one symbol for each expression.")

        # Store each new scalar subexpression in W at the index of its symbol
        handled = set()
        for s, w in zip(vs, ws):
            if W[s] is None:
                W[s] = w
                handled.add(s)
            else:
                assert s in handled  # Result of symmetry!

    # Find symbols of final v from input graph
    vs = V_symbols[-1][0]
    scalar_expression = W[vs]
    return scalar_expression


def _count_nodes_with_unique_post_traversal(expr, e2i=None, skip_terminal_modifiers=False):
    """Yields o for each node o in expr, child before parent.
    Never visits a node twice."""
    if e2i is None:
        e2i = {}

    def getops(e):
        """Get a modifiable list of operands of e, optionally treating modified terminals as a unit."""
        # TODO: Maybe use e._ufl_is_terminal_modifier_
        if e._ufl_is_terminal_ or (skip_terminal_modifiers and is_modified_terminal(e)):
            return []
        else:
            return list(e.ufl_operands)

    stack = [(expr, getops(expr))]
    while stack:
        expr, ops = stack[-1]
        for i, o in enumerate(ops):
            if o is not None and o not in e2i:
                stack.append((o, getops(o)))
                ops[i] = None
                break
        else:
            if not isinstance(expr, (ufl.classes.MultiIndex, ufl.classes.Label)):
                count = len(e2i)
                e2i[expr] = count
            stack.pop()
    return e2i
