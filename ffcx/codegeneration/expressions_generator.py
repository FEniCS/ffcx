# Copyright (C) 2019 Michal Habera
#
# This file is part of FFCX.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import collections
import logging

import ufl
from ffcx.codegeneration.backend import FFCXBackend
from ffcx.codegeneration.C.cnodes import pad_innermost_dim
from ffcx.codegeneration.C.format_lines import format_indented_lines
from ffcx.ir.representationutils import initialize_expression_code

logger = logging.getLogger(__name__)


def generate_expression_code(ir, parameters):

    backend = FFCXBackend(ir, parameters)
    eg = ExpressionGenerator(ir, backend)
    code = initialize_expression_code(ir)
    parts = eg.generate()

    body = format_indented_lines(parts.cs_format(), 1)
    code["tabulate_expression"] = body

    code["original_coefficient_positions"] = format_indented_lines(
        eg.generate_original_coefficient_positions().cs_format(), 1)

    code["points"] = format_indented_lines(eg.generate_points().cs_format(), 1)
    code["value_shape"] = format_indented_lines(eg.generate_value_shape().cs_format(), 1)

    return code


class ExpressionGenerator:
    def __init__(self, ir, backend):

        if len(ir.all_num_points) != 1:
            raise RuntimeError("Only one set of points allowed for expression evaluation")

        self.ir = ir
        self.backend = backend
        self.scope = {}
        self._ufl_names = set()
        self.finalization_blocks = collections.defaultdict(list)
        self.symbol_counters = collections.defaultdict(int)
        self.shared_symbols = {}
        self.num_points = self.ir.all_num_points[0]

    def generate(self):
        L = self.backend.language

        parts = []

        parts += self.generate_element_tables()
        parts += self.generate_unstructured_piecewise_partition()

        all_preparts = []
        all_quadparts = []
        all_postparts = []

        preparts, quadparts, postparts = self.generate_quadrature_loop()
        all_preparts += preparts
        all_quadparts += quadparts
        all_postparts += postparts

        preparts, quadparts, postparts = self.generate_dofblock_partition(quadrature_independent=True)
        all_preparts += preparts
        all_quadparts += quadparts
        all_postparts += postparts

        all_finalizeparts = []

        # Initialize a tensor to zeros
        A_values = [0.0] * ufl.product(self.ir.expression_shape + [self.num_points] + self.ir.tensor_shape)
        all_finalizeparts = self.generate_tensor_value_initialization(A_values)

        # Generate code to add reusable blocks B* to element tensor A
        all_finalizeparts += self.generate_copyout_statements()

        # Collect parts before, during, and after quadrature loops
        parts += all_preparts
        parts += all_quadparts
        parts += all_postparts
        parts += all_finalizeparts

        return L.StatementList(parts)

    def generate_element_tables(self):
        L = self.backend.language
        parts = []

        tables = self.ir.unique_tables

        alignas = self.ir.params["alignas"]
        padlen = self.ir.params["padlen"]
        table_names = sorted(tables)

        for name in table_names:
            table = tables[name]
            decl = L.ArrayDecl(
                "static const ufc_scalar_t", name, table.shape, table, alignas=alignas, padlen=padlen)
            parts += [decl]

        # Add leading comment if there are any tables
        parts = L.commented_code_list(parts, [
            "Precomputed values of basis functions",
            "FE* dimensions: [entities][points][dofs]",
        ])
        return parts

    def generate_quadrature_loop(self):
        """Generate quadrature loop for this num_points."""
        L = self.backend.language

        # Generate unstructured varying partition
        body = self.generate_unstructured_varying_partition()
        body = L.commented_code_list(
            body, "Quadrature loop body setup (num_points={0})".format(self.num_points))

        # Generate dofblock parts, some of this
        # will be placed before or after quadloop
        preparts, quadparts, postparts = \
            self.generate_dofblock_partition()
        body += quadparts

        # Wrap body in loop or scope
        if not body:
            # Could happen for integral with everything zero and optimized away
            quadparts = []
        else:
            iq = self.backend.symbols.quadrature_loop_index()
            quadparts = [L.ForRange(iq, 0, self.num_points, body=body)]

        return preparts, quadparts, postparts

    def generate_unstructured_varying_partition(self):
        L = self.backend.language

        # Get annotated graph of factorisation
        F = self.ir.varying_irs[self.num_points]["factorization"]

        arraysymbol = L.Symbol("sv%d" % self.num_points)
        parts = self.generate_partition(arraysymbol, F, "varying", self.num_points)
        parts = L.commented_code_list(parts, "Unstructured varying computations for num_points=%d" %
                                      (self.num_points, ))
        return parts

    def generate_unstructured_piecewise_partition(self):
        L = self.backend.language

        # Get annotated graph of factorisation
        F = self.ir.piecewise_ir["factorization"]

        arraysymbol = L.Symbol("sp")
        num_points = None
        parts = self.generate_partition(arraysymbol, F, "piecewise", num_points)
        parts = L.commented_code_list(parts, "Unstructured piecewise computations")
        return parts

    def generate_dofblock_partition(self, quadrature_independent=False):
        if quadrature_independent is True:  # NB! None meaning piecewise partition, not custom integral
            block_contributions = self.ir.piecewise_ir["block_contributions"]
        else:
            block_contributions = self.ir.varying_irs[self.num_points]["block_contributions"]

        preparts = []
        quadparts = []
        postparts = []

        blocks = [(blockmap, blockdata)
                  for blockmap, contributions in sorted(block_contributions.items())
                  for blockdata in contributions]

        for blockmap, blockdata in blocks:

            # Define code for block depending on mode
            B, block_preparts, block_quadparts, block_postparts = \
                self.generate_block_parts(self.num_points, blockmap, blockdata, quadrature_independent)

            # Add definitions
            preparts.extend(block_preparts)

            # Add computations
            quadparts.extend(block_quadparts)

            # Add finalization
            postparts.extend(block_postparts)

            # Add A[blockmap] += B[...] to finalization
            self.finalization_blocks[blockmap].append(B)

        return preparts, quadparts, postparts

    def generate_block_parts(self, num_points, blockmap, blockdata, quadrature_independent=False):
        """Generate and return code parts for a given block.

        Returns parts occuring before, inside, and after
        the quadrature loop identified by num_points.

        """
        L = self.backend.language

        # The parts to return
        preparts = []
        quadparts = []
        postparts = []

        alignas = self.ir.params["alignas"]
        padlen = self.ir.params["padlen"]

        block_rank = len(blockmap)
        blockdims = tuple(len(dofmap) for dofmap in blockmap)
        padded_blockdims = pad_innermost_dim(blockdims, padlen)

        ttypes = blockdata.ttypes
        if "zeros" in ttypes:
            raise RuntimeError("Not expecting zero arguments to be left in dofblock generation.")

        iq = self.backend.symbols.quadrature_loop_index()

        arg_indices = tuple(self.backend.symbols.argument_loop_index(i) for i in range(block_rank))

        # Get factor expression
        if blockdata.factor_is_piecewise:
            F = self.ir.piecewise_ir["factorization"]
        else:
            F = self.ir.varying_irs[num_points]["factorization"]

        assert blockdata.block_mode == "full"
        assert not blockdata.transposed, "Not handled yet"

        components = ufl.product(self.ir.expression_shape)

        # Prepend dimensions of dofmap block with free index
        # for quadrature points and expression components
        blockdims = (components, ) + (num_points, ) + blockdims

        B = self.new_temp_symbol("B")
        # Add initialization of this block to parts
        # For all modes, block definition occurs before quadloop
        preparts.append(L.Comment("B[components for block][points][dofs][dofs]"))
        preparts.append(L.ArrayDecl("ufc_scalar_t", B, blockdims, 0, alignas=alignas, padlen=padlen))

        B_indices = tuple([iq] + list(arg_indices))

        # Fetch code to access modified arguments
        # An access to FE table data
        arg_factors = self.get_arg_factors(blockdata, block_rank, num_points, iq, B_indices)

        # Multiply collected factors
        # A list of computations of Bs, for each component of the factor expression
        # Add result to block inside quadloop
        body = []

        for fi_ci in blockdata.factor_indices_comp_indices:
            f = self.get_var(num_points, F.nodes[fi_ci[0]]["expression"])
            Brhs = L.float_product([f] + arg_factors)
            body.append(L.AssignAdd(B[(fi_ci[1],) + B_indices], Brhs))

        for i in reversed(range(block_rank)):
            body = L.ForRange(
                B_indices[i + 1], 0, padded_blockdims[i], body=body)
        quadparts += [body]

        # Define rhs expression for A[it][iq][blockmap[arg_indices]] += A_rhs
        # This is used outside quadloop, in finalization copyout statements
        iec = self.backend.symbols.expr_component_index()
        A_rhs = B[(iec,) + B_indices]

        return A_rhs, preparts, quadparts, postparts

    def generate_tensor_value_initialization(self, A_values):
        parts = []

        L = self.backend.language
        A = self.backend.symbols.element_tensor()
        A_size = len(A_values)

        z = L.LiteralFloat(0.0)

        k = L.Symbol("k")

        # Zero everything first
        parts += [L.ForRange(k, 0, A_size, index_type="int", body=L.Assign(A[k], 0.0))]

        # Generate A[i] = A_values[i] skipping zeros
        for i in range(A_size):
            if not (A_values[i] == 0.0 or A_values[i] == z):
                parts += [L.Assign(A[i], A_values[i])]

        return parts

    def generate_copyout_statements(self):
        L = self.backend.language
        parts = []

        # Get symbol, dimensions, and loop index symbols for A
        A_shape = self.ir.tensor_shape
        A_rank = len(A_shape)

        Asym = self.backend.symbols.element_tensor()

        num_expr_components = ufl.product(self.ir.expression_shape)
        A = L.FlattenedArray(Asym, dims=[num_expr_components] + [self.num_points] + A_shape)

        indices = [self.backend.symbols.argument_loop_index(i) for i in range(A_rank)]

        dofmap_parts = []
        dofmaps = {}
        for blockmap, contributions in sorted(self.finalization_blocks.items()):

            # Define mapping from B indices to A indices
            A_indices = []
            for i in range(A_rank):
                dofmap = blockmap[i]
                begin = dofmap[0]
                end = dofmap[-1] + 1
                if len(dofmap) == end - begin:
                    # Dense insertion, offset B index to index A
                    j = indices[i] + begin
                else:
                    # Sparse insertion, map B index through dofmap
                    DM = dofmaps.get(dofmap)
                    if DM is None:
                        DM = L.Symbol("DM%d" % len(dofmaps))
                        dofmaps[dofmap] = DM
                        dofmap_parts.append(
                            L.ArrayDecl("static const int", DM, len(dofmap), dofmap))
                    j = DM[indices[i]]
                A_indices.append(j)

            # Sum up all blocks contributing to this blockmap
            term = L.Sum([B_rhs for B_rhs in contributions])

            # Add components of all B's to A component in loop nest
            iq = self.backend.symbols.quadrature_loop_index()
            iec = self.backend.symbols.expr_component_index()
            body = L.AssignAdd(A[tuple([iec] + [iq] + A_indices)], term)

            for i in reversed(range(A_rank)):
                index_range = len(blockmap[i])
                body = L.ForRange(indices[i], 0, index_range, body=body)

            body = L.ForRange(iq, 0, self.num_points, body=body)

            # TODO: Here we assume that block contributes to each
            #       component of the expression. In some cases produces
            #       suboptimal copyout statements (adding zeros to tensor A)
            body = L.ForRange(iec, 0, num_expr_components, body=body)

            # Add this block to parts
            parts.append(body)

        # Place static dofmap tables first
        parts = dofmap_parts + parts

        parts = L.commented_code_list(parts, "Copyout blocks B into A, respecting dofmaps")
        return parts

    def get_arg_factors(self, blockdata, block_rank, num_points, iq, indices):
        L = self.backend.language

        arg_factors = []
        for i in range(block_rank):
            mad = blockdata.ma_data[i]
            td = mad.tabledata
            mt = self.ir.piecewise_ir["modified_arguments"][mad.ma_index]

            table = self.backend.symbols.element_table(td, self.ir.entitytype, mt.restriction)

            assert td.ttype != "zeros"

            if td.ttype == "ones":
                arg_factor = L.LiteralFloat(1.0)
            else:
                # Assuming B sparsity follows element table sparsity
                arg_factor = table[indices[i + 1]]
            arg_factors.append(arg_factor)
        return arg_factors

    def new_temp_symbol(self, basename):
        """Create a new code symbol named basename + running counter."""
        L = self.backend.language
        name = "%s%d" % (basename, self.symbol_counters[basename])
        self.symbol_counters[basename] += 1
        return L.Symbol(name)

    def get_var(self, num_points, v):
        if v._ufl_is_literal_:
            return self.backend.ufl_to_language.get(v)
        f = self.scope.get(v)
        return f

    def generate_partition(self, symbol, F, mode, num_points):
        L = self.backend.language

        definitions = []
        intermediates = []

        for i, attr in F.nodes.items():
            if attr['status'] != mode:
                continue
            v = attr['expression']
            mt = attr.get('mt')

            if v._ufl_is_literal_:
                vaccess = self.backend.ufl_to_language.get(v)
            elif mt is not None:
                # All finite element based terminals have table data, as well
                # as some, but not all, of the symbolic geometric terminals
                tabledata = attr.get('tr')

                # Backend specific modified terminal translation
                vaccess = self.backend.access.get(mt.terminal, mt, tabledata, num_points)
                vdef = self.backend.definitions.get(mt.terminal, mt, tabledata, num_points, vaccess)

                # Store definitions of terminals in list
                assert isinstance(vdef, list)
                definitions.extend(vdef)
            else:
                # Get previously visited operands
                vops = [self.get_var(num_points, op) for op in v.ufl_operands]

                # get parent operand
                pid = F.in_edges[i][0] if F.in_edges[i] else -1
                if pid and pid > i:
                    parent_exp = F.nodes.get(pid)['expression']
                else:
                    parent_exp = None

                # Mapping UFL operator to target language
                self._ufl_names.add(v._ufl_handler_name_)
                vexpr = self.backend.ufl_to_language.get(v, *vops)

                # Create a new intermediate for each subexpression
                # except boolean conditions and its childs
                if isinstance(parent_exp, ufl.classes.Condition):
                    # Skip intermediates for 'x' and 'y' in x<y
                    # Avoid the creation of complex valued intermediates
                    vaccess = vexpr
                elif isinstance(v, ufl.classes.Condition):
                    # Inline the conditions x < y, condition values
                    # This removes the need to handle boolean intermediate variables.
                    # With tensor-valued conditionals it may not be optimal but we
                    # let the compiler take responsibility for optimizing those cases.
                    vaccess = vexpr
                elif any(op._ufl_is_literal_ for op in v.ufl_operands):
                    # Skip intermediates for e.g. -2.0*x,
                    # resulting in lines like z = y + -2.0*x
                    vaccess = vexpr
                else:
                    # Record assignment of vexpr to intermediate variable
                    j = len(intermediates)
                    if self.ir.params["use_symbol_array"]:
                        vaccess = symbol[j]
                        intermediates.append(L.Assign(vaccess, vexpr))
                    else:
                        vaccess = L.Symbol("%s_%d" % (symbol.name, j))
                        intermediates.append(L.VariableDecl("const ufc_scalar_t", vaccess, vexpr))

            # Store access node for future reference
            self.scope[v] = vaccess

        # Join terminal computation, array of intermediate expressions,
        # and intermediate computations
        parts = []
        if definitions:
            parts += definitions
        if intermediates:
            if self.ir.params["use_symbol_array"]:
                alignas = self.ir.params["alignas"]
                parts += [L.ArrayDecl("ufc_scalar_t", symbol, len(intermediates), alignas=alignas)]
            parts += intermediates
        return parts

    def generate_original_coefficient_positions(self):
        L = self.backend.language
        num_coeffs = len(self.ir.original_coefficient_positions)
        orig_pos = L.Symbol("original_coefficient_positions")
        if num_coeffs > 0:
            parts = [L.ArrayDecl("static const int", orig_pos,
                                 values=self.ir.original_coefficient_positions,
                                 sizes=(num_coeffs, ))]
            parts += [L.Assign("expression->original_coefficient_positions", orig_pos)]
        else:
            parts = []
        return L.StatementList(parts)

    def generate_points(self):
        L = self.backend.language
        parts = L.ArrayDecl("static const double", "points", values=self.ir.points,
                            sizes=self.ir.points.shape)
        return parts

    def generate_value_shape(self):
        L = self.backend.language
        parts = L.ArrayDecl("static const int", "value_shape", values=self.ir.expression_shape,
                            sizes=len(self.ir.expression_shape))
        return parts
