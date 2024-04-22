# Copyright (C) 2011-2017 Martin Sandve Alnæs
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Algorithms for value numbering within computational graphs."""

import logging

import ufl
from ufl.pullback import SymmetricPullback

from ffcx.ir.analysis.indexing import (
    map_component_tensor_arg_components,
    map_indexed_arg_components,
)
from ffcx.ir.analysis.modified_terminals import analyse_modified_terminal

logger = logging.getLogger("ffcx")


class ValueNumberer:
    """Maps scalar components to unique values.

    An algorithm to map the scalar components of an expression node to unique value numbers,
    with fallthrough for types that can be mapped to the value numbers
    of their operands.
    """

    def __init__(self, G):
        """Initialise."""
        self.symbol_count = 0
        self.G = G
        self.V_symbols = []
        self.call_lookup = {
            ufl.classes.Expr: self.expr,
            ufl.classes.Argument: self.form_argument,
            ufl.classes.Coefficient: self.form_argument,
            ufl.classes.Grad: self._modified_terminal,
            ufl.classes.ReferenceGrad: self._modified_terminal,
            ufl.classes.FacetAvg: self._modified_terminal,
            ufl.classes.CellAvg: self._modified_terminal,
            ufl.classes.Restricted: self._modified_terminal,
            ufl.classes.ReferenceValue: self._modified_terminal,
            ufl.classes.Indexed: self.indexed,
            ufl.classes.ComponentTensor: self.component_tensor,
            ufl.classes.ListTensor: self.list_tensor,
            ufl.classes.Variable: self.variable,
        }

    def new_symbols(self, n):
        """Generate new symbols with a running counter."""
        begin = self.symbol_count
        end = begin + n
        self.symbol_count = end
        return list(range(begin, end))

    def new_symbol(self):
        """Generate new symbol with a running counter."""
        begin = self.symbol_count
        self.symbol_count += 1
        return begin

    def get_node_symbols(self, expr):
        """Get node symbols."""
        idx = [i for i, v in self.G.nodes.items() if v["expression"] == expr][0]
        return self.V_symbols[idx]

    def compute_symbols(self):
        """Compute symbols."""
        for i, v in self.G.nodes.items():
            expr = v["expression"]
            symbol = None
            # First look for exact type match
            f = self.call_lookup.get(type(expr), False)
            if f:
                symbol = f(expr)
            else:
                # Look for parent class types instead
                for k in self.call_lookup.keys():
                    if isinstance(expr, k):
                        symbol = self.call_lookup[k](expr)
                        break

            if symbol is None:
                # Nothing found
                raise RuntimeError("Not expecting type %s here." % type(expr))

            self.V_symbols.append(symbol)

        return self.V_symbols

    def expr(self, v):
        """Create new symbols for expressions that represent new values."""
        n = ufl.product(v.ufl_shape + v.ufl_index_dimensions)
        return self.new_symbols(n)

    def form_argument(self, v):
        """Create new symbols for expressions that represent new values."""
        e = v.ufl_function_space().ufl_element()

        if isinstance(e.pullback, SymmetricPullback):
            # Build symbols with symmetric components skipped
            symbols = []
            mapped_symbols = {}
            for c in ufl.permutation.compute_indices(v.ufl_shape):
                # Build mapped component mc with symmetries from element considered
                mc = min(i for i, j in e.pullback._symmetry.items() if j == e.pullback._symmetry[c])

                # Get existing symbol or create new and store with mapped component mc as key
                s = mapped_symbols.get(mc)
                if s is None:
                    s = self.new_symbol()
                    mapped_symbols[mc] = s
                symbols.append(s)
        else:
            n = ufl.product(v.ufl_shape + v.ufl_index_dimensions)
            symbols = self.new_symbols(n)

        return symbols

    # Handle modified terminals with element symmetries and second derivative symmetries!
    # terminals are implemented separately, or maybe they don't need to be?

    def _modified_terminal(self, v):
        """Handle modified terminal.

        Modifiers:
            terminal: the underlying Terminal object
            global_derivatives: tuple of ints, each meaning derivative
                in that global direction
            local_derivatives: tuple of ints, each meaning derivative in
                that local direction
            reference_value: bool, whether this is represented in
                reference frame
            averaged: None, 'facet' or 'cell'
            restriction: None, '+' or '-'
            component: tuple of ints, the global component of the Terminal
                flat_component: single int, flattened local component of the
            Terminal, considering symmetry
        """
        # (1) mt.terminal.ufl_shape defines a core indexing space UNLESS mt.reference_value,
        #     in which case the reference value shape of the element must be used.
        # (2) mt.terminal.ufl_element().symmetry() defines core symmetries
        # (3) averaging and restrictions define distinct symbols, no additional symmetries
        # (4) two or more grad/reference_grad defines distinct symbols with additional symmetries

        # v is not necessary scalar here, indexing in (0,...,0) picks the first scalar component
        # to analyse, which should be sufficient to get the base shape and derivatives
        if v.ufl_shape:
            mt = analyse_modified_terminal(v[(0,) * len(v.ufl_shape)])
        else:
            mt = analyse_modified_terminal(v)

        # Get derivatives
        num_ld = len(mt.local_derivatives)
        num_gd = len(mt.global_derivatives)
        assert not (num_ld and num_gd)
        if num_ld:
            domain = ufl.domain.extract_unique_domain(mt.terminal)
            tdim = domain.topological_dimension()
            d_components = ufl.permutation.compute_indices((tdim,) * num_ld)
        elif num_gd:
            domain = ufl.domain.extract_unique_domiain(mt.terminal)
            gdim = domain.geometric_dimension()
            d_components = ufl.permutation.compute_indices((gdim,) * num_gd)
        else:
            d_components = [()]

        # Get base shape without the derivative axes
        base_components = ufl.permutation.compute_indices(mt.base_shape)

        # Build symbols with symmetric components and derivatives skipped
        symbols = []
        mapped_symbols = {}
        for bc in base_components:
            for dc in d_components:
                # Build mapped component mc with symmetries from element
                # and derivatives combined
                mbc = mt.base_symmetry.get(bc, bc)
                mdc = tuple(sorted(dc))
                mc = mbc + mdc

                # Get existing symbol or create new and store with
                # mapped component mc as key
                s = mapped_symbols.get(mc)
                if s is None:
                    s = self.new_symbol()
                    mapped_symbols[mc] = s
                symbols.append(s)

        # Consistency check before returning symbols
        assert not v.ufl_free_indices
        if ufl.product(v.ufl_shape) != len(symbols):
            raise RuntimeError("Internal error in value numbering.")
        return symbols

    def indexed(self, Aii):
        """Return indexed value.

        This is implemented as a fall-through operation.
        """
        # Reuse symbols of arg A for Aii
        A = Aii.ufl_operands[0]

        # Get symbols of argument A
        A_symbols = self.get_node_symbols(A)

        # Map A_symbols to Aii_symbols
        d = map_indexed_arg_components(Aii)
        symbols = [A_symbols[k] for k in d]
        return symbols

    def component_tensor(self, A):
        """Component tensor."""
        # Reuse symbols of arg Aii for A
        Aii = A.ufl_operands[0]

        # Get symbols of argument Aii
        Aii_symbols = self.get_node_symbols(Aii)

        # Map A_symbols to Aii_symbols
        d = map_component_tensor_arg_components(A)
        symbols = [Aii_symbols[k] for k in d]
        return symbols

    def list_tensor(self, v):
        """List tensor."""
        symbols = []
        for row in v.ufl_operands:
            symbols.extend(self.get_node_symbols(row))
        return symbols

    def variable(self, v):
        """Direct reuse of all symbols."""
        return self.get_node_symbols(v.ufl_operands[0])
