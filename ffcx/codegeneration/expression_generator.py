# Copyright (C) 2019 Michal Habera
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Expression generator."""

import collections
import logging
from itertools import product
from typing import Any

import ufl

import ffcx.codegeneration.lnodes as L
from ffcx.codegeneration import geometry
from ffcx.codegeneration.backend import FFCXBackend
from ffcx.codegeneration.lnodes import LNode
from ffcx.ir.representation import ExpressionIR

logger = logging.getLogger("ffcx")


class ExpressionGenerator:
    """Expression generator."""

    def __init__(self, ir: ExpressionIR, backend: FFCXBackend):
        """Initialise."""
        if len(list(ir.integrand.keys())) != 1:
            raise RuntimeError("Only one set of points allowed for expression evaluation")

        self.ir = ir
        self.backend = backend
        self.scope: dict[Any, LNode] = {}
        self._ufl_names: set[Any] = set()
        self.symbol_counters: collections.defaultdict[Any, int] = collections.defaultdict(int)
        self.shared_symbols: dict[Any, Any] = {}
        self.quadrature_rule = list(self.ir.integrand.keys())[0]

    def generate(self):
        """Generate."""
        parts = []
        parts += self.generate_element_tables()

        # Generate the tables of geometry data that are needed
        parts += self.generate_geometry_tables()
        parts += self.generate_piecewise_partition()

        all_preparts = []
        all_quadparts = []

        preparts, quadparts = self.generate_quadrature_loop()
        all_preparts += preparts
        all_quadparts += quadparts

        # Collect parts before, during, and after quadrature loops
        parts += all_preparts
        parts += all_quadparts

        return L.StatementList(parts)

    def generate_geometry_tables(self):
        """Generate static tables of geometry data."""
        # Currently we only support circumradius
        ufl_geometry = {
            ufl.geometry.ReferenceCellVolume: "reference_cell_volume",
            ufl.geometry.ReferenceNormal: "reference_facet_normals",
        }

        cells: dict[Any, set[Any]] = {t: set() for t in ufl_geometry.keys()}  # type: ignore
        for integrand in self.ir.integrand.values():
            for attr in integrand["factorization"].nodes.values():
                mt = attr.get("mt")
                if mt is not None:
                    t = type(mt.terminal)
                    if self.ir.entitytype == "cell" and issubclass(
                        t, ufl.geometry.GeometricFacetQuantity
                    ):
                        raise RuntimeError(f"Expressions for cells do not support {t}.")
                    if t in ufl_geometry:
                        cells[t].add(
                            ufl.domain.extract_unique_domain(mt.terminal).ufl_cell().cellname()
                        )

        parts = []
        for i, cell_list in cells.items():
            for c in cell_list:
                parts.append(geometry.write_table(ufl_geometry[i], c))

        return parts

    def generate_element_tables(self):
        """Generate tables of FE basis evaluated at specified points."""
        parts = []

        tables = self.ir.unique_tables
        table_names = sorted(tables)

        for name in table_names:
            table = tables[name]
            symbol = L.Symbol(name, dtype=L.DataType.REAL)
            self.backend.symbols.element_tables[name] = symbol
            decl = L.ArrayDecl(symbol, sizes=table.shape, values=table, const=True)
            parts += [decl]

        # Add leading comment if there are any tables
        parts = L.commented_code_list(
            parts,
            [
                "Precomputed values of basis functions",
                "FE* dimensions: [entities][points][dofs]",
            ],
        )
        return parts

    def generate_quadrature_loop(self):
        """Generate quadrature loop for this quadrature rule.

        In the context of expressions quadrature loop is not accumulated.
        """
        # Generate varying partition
        body = self.generate_varying_partition()
        body = L.commented_code_list(
            body, f"Points loop body setup quadrature loop {self.quadrature_rule.id()}"
        )

        # Generate dofblock parts, some of this
        # will be placed before or after quadloop
        preparts, quadparts = self.generate_dofblock_partition()
        body += quadparts

        # Wrap body in loop or scope
        if not body:
            # Could happen for integral with everything zero and optimized away
            quadparts = []
        else:
            iq = self.backend.symbols.quadrature_loop_index
            num_points = self.quadrature_rule.points.shape[0]
            quadparts = [L.ForRange(iq, 0, num_points, body=body)]
        return preparts, quadparts

    def generate_varying_partition(self):
        """Generate factors of blocks which are not cellwise constant."""
        # Get annotated graph of factorisation
        F = self.ir.integrand[self.quadrature_rule]["factorization"]

        arraysymbol = L.Symbol(f"sv_{self.quadrature_rule.id()}", dtype=L.DataType.SCALAR)
        parts = self.generate_partition(arraysymbol, F, "varying")
        parts = L.commented_code_list(
            parts,
            f"Unstructured varying computations for quadrature rule {self.quadrature_rule.id()}",
        )
        return parts

    def generate_piecewise_partition(self):
        """Generate factors of blocks which are constant.

        I.e. do not depend on quadrature points).
        """
        # Get annotated graph of factorisation
        F = self.ir.integrand[self.quadrature_rule]["factorization"]

        arraysymbol = L.Symbol("sp", dtype=L.DataType.SCALAR)
        parts = self.generate_partition(arraysymbol, F, "piecewise")
        parts = L.commented_code_list(parts, "Unstructured piecewise computations")
        return parts

    def generate_dofblock_partition(self):
        """Generate assignments of blocks multiplied with their factors into final tensor A."""
        block_contributions = self.ir.integrand[self.quadrature_rule]["block_contributions"]

        preparts = []
        quadparts = []

        blocks = [
            (blockmap, blockdata)
            for blockmap, contributions in sorted(block_contributions.items())
            for blockdata in contributions
        ]

        for blockmap, blockdata in blocks:
            # Define code for block depending on mode
            block_preparts, block_quadparts = self.generate_block_parts(blockmap, blockdata)

            # Add definitions
            preparts.extend(block_preparts)

            # Add computations
            quadparts.extend(block_quadparts)

        return preparts, quadparts

    def generate_block_parts(self, blockmap, blockdata):
        """Generate and return code parts for a given block."""
        # The parts to return
        preparts = []
        quadparts = []

        block_rank = len(blockmap)
        blockdims = tuple(len(dofmap) for dofmap in blockmap)

        ttypes = blockdata.ttypes
        if "zeros" in ttypes:
            raise RuntimeError("Not expecting zero arguments to be left in dofblock generation.")

        arg_indices = tuple(self.backend.symbols.argument_loop_index(i) for i in range(block_rank))

        F = self.ir.integrand[self.quadrature_rule]["factorization"]

        assert not blockdata.transposed, "Not handled yet"
        components = ufl.product(self.ir.expression_shape)

        num_points = self.quadrature_rule.points.shape[0]
        A_shape = [num_points, components] + self.ir.tensor_shape
        A = self.backend.symbols.element_tensor
        iq = self.backend.symbols.quadrature_loop_index

        # Check if DOFs in dofrange are equally spaced.
        expand_loop = False
        for bm in blockmap:
            for a, b in zip(bm[1:-1], bm[2:]):
                if b - a != bm[1] - bm[0]:
                    expand_loop = True
                    break
            else:
                continue
            break

        if expand_loop:
            # If DOFs in dofrange are not equally spaced, then expand out the for loop
            for A_indices, B_indices in zip(
                product(*blockmap), product(*[range(len(b)) for b in blockmap])
            ):
                B_indices = tuple([iq] + list(B_indices))
                A_indices = tuple([iq] + A_indices)
                for fi_ci in blockdata.factor_indices_comp_indices:
                    f = self.get_var(F.nodes[fi_ci[0]]["expression"])
                    arg_factors = self.get_arg_factors(blockdata, block_rank, B_indices)
                    Brhs = L.float_product([f] + arg_factors)
                    multi_index = L.MultiIndex([A_indices[0], fi_ci[1]] + A_indices[1:], A_shape)
                    quadparts.append(L.AssignAdd(A[multi_index], Brhs))
        else:
            # Prepend dimensions of dofmap block with free index
            # for quadrature points and expression components
            B_indices = tuple([iq] + list(arg_indices))

            # Fetch code to access modified arguments
            # An access to FE table data
            arg_factors = self.get_arg_factors(blockdata, block_rank, B_indices)

            # TODO: handle non-contiguous dof ranges

            A_indices = []
            for bm, index in zip(blockmap, arg_indices):
                # TODO: switch order here? (optionally)
                offset = bm[0]
                if len(bm) == 1:
                    A_indices.append(index + offset)
                else:
                    block_size = bm[1] - bm[0]
                    A_indices.append(block_size * index + offset)
            A_indices = tuple([iq] + A_indices)

            # Multiply collected factors
            # For each component of the factor expression
            # add result inside quadloop
            body = []

            for fi_ci in blockdata.factor_indices_comp_indices:
                f = self.get_var(F.nodes[fi_ci[0]]["expression"])
                Brhs = L.float_product([f] + arg_factors)
                indices = [A_indices[0], fi_ci[1]] + list(A_indices[1:])
                multi_index = L.MultiIndex(indices, A_shape)
                body.append(L.AssignAdd(A[multi_index], Brhs))

            for i in reversed(range(block_rank)):
                body = L.ForRange(B_indices[i + 1], 0, blockdims[i], body=body)
            quadparts += [body]

        return preparts, quadparts

    def get_arg_factors(self, blockdata, block_rank, indices):
        """Get argument factors (i.e. blocks).

        Args:
            blockdata: block data
            block_rank: block rank
            indices: Indices used to index element tables
        """
        arg_factors = []
        for i in range(block_rank):
            mad = blockdata.ma_data[i]
            td = mad.tabledata
            mt = self.ir.integrand[self.quadrature_rule]["modified_arguments"][mad.ma_index]

            table = self.backend.symbols.element_table(td, self.ir.entitytype, mt.restriction)

            assert td.ttype != "zeros"

            if td.ttype == "ones":
                arg_factor = L.LiteralFloat(1.0)
            else:
                arg_factor = table[indices[i + 1]]
            arg_factors.append(arg_factor)
        return arg_factors

    def new_temp_symbol(self, basename):
        """Create a new code symbol named basename + running counter."""
        name = "%s%d" % (basename, self.symbol_counters[basename])
        self.symbol_counters[basename] += 1
        return L.Symbol(name, dtype=L.DataType.SCALAR)

    def get_var(self, v):
        """Get a variable."""
        if v._ufl_is_literal_:
            return L.ufl_to_lnodes(v)
        f = self.scope.get(v)
        return f

    def generate_partition(self, symbol, F, mode):
        """Generate computations of factors of blocks."""
        definitions = []
        intermediates = []

        for _, attr in F.nodes.items():
            if attr["status"] != mode:
                continue
            v = attr["expression"]
            mt = attr.get("mt")

            if v._ufl_is_literal_:
                vaccess = L.ufl_to_lnodes(v)
            elif mt is not None:
                # All finite element based terminals have table data, as well
                # as some, but not all, of the symbolic geometric terminals
                tabledata = attr.get("tr")

                # Backend specific modified terminal translation
                vaccess = self.backend.access.get(mt, tabledata, 0)
                vdef = self.backend.definitions.get(mt, tabledata, 0, vaccess)

                if vdef:
                    assert isinstance(vdef, L.Section)
                    vdef = vdef.declarations + vdef.statements

                # Store definitions of terminals in list
                assert isinstance(vdef, list)
                definitions.extend(vdef)
            else:
                # Get previously visited operands
                vops = [self.get_var(op) for op in v.ufl_operands]

                # Mapping UFL operator to target language
                self._ufl_names.add(v._ufl_handler_name_)
                vexpr = L.ufl_to_lnodes(v, *vops)

                is_cond = isinstance(v, ufl.classes.Condition)
                dtype = L.DataType.BOOL if is_cond else L.DataType.SCALAR

                j = len(intermediates)
                vaccess = L.Symbol(f"{symbol.name}_{j}", dtype=dtype)
                intermediates.append(L.VariableDecl(vaccess, vexpr))

            # Store access node for future reference
            self.scope[v] = vaccess

        # Join terminal computation, array of intermediate expressions,
        # and intermediate computations
        parts = []

        parts += definitions
        parts += intermediates

        return parts
