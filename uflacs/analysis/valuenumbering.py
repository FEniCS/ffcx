
from six.moves import xrange as range
from ufl.common import product
from ufl.permutation import compute_indices
from ufl.algorithms import MultiFunction
from ufl.classes import FormArgument
from uflacs.analysis.indexing import map_indexed_arg_components, map_component_tensor_arg_components
from uflacs.analysis.modified_terminals import analyse_modified_terminal


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
        symmetry = v.element().symmetry()

        if False and symmetry:
            # FIXME: Ignoring symmetries for now, handle by creating only
            # some new symbols and mapping the rest using the symmetry map.
            actual_components = sorted(set(symmetry.values()))
            m = len(actual_components)
            actual_symbols = self.new_symbols(m)
            symbols = mapping_of_actual_symbols_to_all_components(actual_symbols, symmetry)  # Need to implement this

        else:
            n = self.V_sizes[i]
            symbols = self.new_symbols(n)

        return symbols

    def _modified_terminal(self, v, i):
        """Modifiers:
        terminal           - the underlying Terminal object
        global_derivatives - tuple of ints, each meaning derivative in that global direction
        local_derivatives  - tuple of ints, each meaning derivative in that local direction
        averaged           - None, 'facet' or 'cell'
        restriction        - None, '+' or '-'
        component          - tuple of ints, the global component of the Terminal
        flat_component     - single int, flattened local component of the Terminal, considering symmetry
        """

        # (1) mt.terminal.shape() defines a core indexing space
        # (2) mt.terminal.element().symmetry() defines core symmetries
        # (3) averaging and restrictions define distinct symbols, no additional symmetries
        # (4) two or more grad/reference_grad defines distinct symbols with additional symmetries

        # FIXME: Need modified version of amt(), v is probably not scalar here. This hack works for now.
        if v.shape():
            mt = analyse_modified_terminal(v[(0,) * len(v.shape())])
        else:
            mt = analyse_modified_terminal(v)

        cell = mt.terminal.cell()
        tdim = cell.topological_dimension()
        gdim = cell.geometric_dimension()

        num_ld = len(mt.local_derivatives)
        num_gd = len(mt.global_derivatives)

        base_components = compute_indices(mt.terminal.shape())
        assert not (num_ld and num_gd)
        if num_ld:
            # d_components = compute_permutations(num_ld, tdim)
            d_components = compute_indices((tdim,) * num_ld)
        elif num_gd:
            # d_components = compute_permutations(num_gd, gdim)
            d_components = compute_indices((gdim,) * num_gd)
        else:
            d_components = [()]

        if isinstance(mt, FormArgument):
            symmetry = mt.element().symmetry()
        else:
            symmetry = {}

        symbols = []
        mapped_symbols = {}
        for bc in base_components:
            for dc in d_components:
                # Build mapped component with symmetries from element and derivatives combined
                mbc = symmetry.get(bc, bc)
                mdc = tuple(sorted(dc))
                c = bc + dc
                mc = mbc + mdc

                # Get existing symbol or create new and store with mapped component mc as key
                s = mapped_symbols.get(mc)
                if s is None:
                    s = self.new_symbol()
                    mapped_symbols[mc] = s
                symbols.append(s)

        assert not v.free_indices()
        if not product(v.shape()) == len(symbols):
            from ufl.algorithms import tree_format
            print()
            print(num_ld)
            print(num_gd)
            print(d_components)
            print(v.shape())
            print(len(symbols))
            print(tree_format(v))
            print()

        return symbols

    # Handle modified terminals with element symmetries and second derivative symmetries!
    # terminals are implemented separately, or maybe they don't need to be?
    grad = _modified_terminal
    reference_grad = _modified_terminal
    facet_avg = _modified_terminal
    cell_avg = _modified_terminal
    restricted = _modified_terminal
    reference_value = _modified_terminal  # Not yet implemented in UFL
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
            symbols.extend(rowsymb)  # FIXME: Test that this produces the right transposition
        return symbols

    def transposed(self, AT, i):
        A, = AT.ufl_operands

        assert not A.free_indices(), "Assuming no free indices in transposed (for now), report as bug if needed."  # FIXME
        r, c = A.shape()

        A_symbols = self.get_node_symbols(A)
        assert len(A_symbols) == r * c

        # AT[j,i] = A[i,j]
        # sh(A) = (r,c)
        # sh(AT) = (c,r)
        # AT[j*r+i] = A[i*c+j]
        symbols = [None] * (r * c)
        for j in range(c):
            for i in range(r):
                symbols[j * r + i] = A_symbols[i * c + j]
        return symbols

    def variable(self, v, i):
        "Direct reuse of all symbols."
        return self.get_node_symbols(v.ufl_operands[0])
