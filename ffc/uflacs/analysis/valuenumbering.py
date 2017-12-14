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

"""Algorithms for value numbering within computational graphs."""

from ufl import product
from ufl.permutation import compute_indices
from ufl.corealg.multifunction import MultiFunction
from ufl.classes import FormArgument

from ffc.log import error

from ffc.uflacs.analysis.indexing import map_indexed_arg_components, map_component_tensor_arg_components
from ffc.uflacs.analysis.modified_terminals import analyse_modified_terminal


class ValueNumberer(MultiFunction):

    """An algorithm to map the scalar components of an expression node to unique value numbers,
    with fallthrough for types that can be mapped to the value numbers of their operands."""

    def __init__(self, e2i, V_sizes, V_symbols):
        MultiFunction.__init__(self)
        self.symbol_count = 0
        self.e2i = e2i
        self.V_sizes = V_sizes
        self.V_symbols = V_symbols

    def new_symbols(self, n):
        "Generator for new symbols with a running counter."
        begin = self.symbol_count
        end = begin + n
        self.symbol_count = end
        return list(range(begin, end))

    def new_symbol(self):
        "Generator for new symbols with a running counter."
        begin = self.symbol_count
        self.symbol_count += 1
        return begin

    def get_node_symbols(self, expr):
        return self.V_symbols[self.e2i[expr]]

    def expr(self, v, i):
        "Create new symbols for expressions that represent new values."
        n = self.V_sizes[i]
        return self.new_symbols(n)

    def form_argument(self, v, i):
        "Create new symbols for expressions that represent new values."
        symmetry = v.ufl_element().symmetry()

        if symmetry:
            # Build symbols with symmetric components skipped
            symbols = []
            mapped_symbols = {}
            for c in compute_indices(v.ufl_shape):
                # Build mapped component mc with symmetries from element considered
                mc = symmetry.get(c, c)

                # Get existing symbol or create new and store with mapped component mc as key
                s = mapped_symbols.get(mc)
                if s is None:
                    s = self.new_symbol()
                    mapped_symbols[mc] = s
                symbols.append(s)
        else:
            n = self.V_sizes[i]
            symbols = self.new_symbols(n)

        return symbols

    def _modified_terminal(self, v, i):
        """Modifiers:
        terminal           - the underlying Terminal object
        global_derivatives - tuple of ints, each meaning derivative in that global direction
        local_derivatives  - tuple of ints, each meaning derivative in that local direction
        reference_value    - bool, whether this is represented in reference frame
        averaged           - None, 'facet' or 'cell'
        restriction        - None, '+' or '-'
        component          - tuple of ints, the global component of the Terminal
        flat_component     - single int, flattened local component of the Terminal, considering symmetry
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
            domain = mt.terminal.ufl_domain()
            tdim = domain.topological_dimension()
            d_components = compute_indices((tdim,) * num_ld)
        elif num_gd:
            domain = mt.terminal.ufl_domain()
            gdim = domain.geometric_dimension()
            d_components = compute_indices((gdim,) * num_gd)
        else:
            d_components = [()]

        # Get base shape without the derivative axes
        base_components = compute_indices(mt.base_shape)

        # Build symbols with symmetric components and derivatives skipped
        symbols = []
        mapped_symbols = {}
        for bc in base_components:
            for dc in d_components:
                # Build mapped component mc with symmetries from element and derivatives combined
                mbc = mt.base_symmetry.get(bc, bc)
                mdc = tuple(sorted(dc))
                mc = mbc + mdc

                # Get existing symbol or create new and store with mapped component mc as key
                s = mapped_symbols.get(mc)
                if s is None:
                    s = self.new_symbol()
                    mapped_symbols[mc] = s
                symbols.append(s)

        # Consistency check before returning symbols
        assert not v.ufl_free_indices
        if product(v.ufl_shape) != len(symbols):
            error("Internal error in value numbering.")
        return symbols

    # Handle modified terminals with element symmetries and second derivative symmetries!
    # terminals are implemented separately, or maybe they don't need to be?
    grad = _modified_terminal
    reference_grad = _modified_terminal
    facet_avg = _modified_terminal
    cell_avg = _modified_terminal
    restricted = _modified_terminal
    reference_value = _modified_terminal
    # indexed is implemented as a fall-through operation

    def indexed(self, Aii, i):
        # Reuse symbols of arg A for Aii
        A = Aii.ufl_operands[0]

        # Get symbols of argument A
        A_symbols = self.get_node_symbols(A)

        # Map A_symbols to Aii_symbols
        d = map_indexed_arg_components(Aii)
        symbols = [A_symbols[k] for k in d]
        return symbols

    def component_tensor(self, A, i):
        # Reuse symbols of arg Aii for A
        Aii = A.ufl_operands[0]

        # Get symbols of argument Aii
        Aii_symbols = self.get_node_symbols(Aii)

        # Map A_symbols to Aii_symbols
        d = map_component_tensor_arg_components(A)
        symbols = [Aii_symbols[k] for k in d]
        return symbols

    def list_tensor(self, v, i):
        row_symbols = [self.get_node_symbols(row) for row in v.ufl_operands]
        symbols = []
        for rowsymb in row_symbols:
            symbols.extend(rowsymb)
        return symbols

    def variable(self, v, i):
        "Direct reuse of all symbols."
        return self.get_node_symbols(v.ufl_operands[0])
