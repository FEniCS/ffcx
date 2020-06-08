# Copyright (C) 2019 Michal Habera
#
# This file is part of FFCX.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
import collections
import logging

import ufl
from ffcx.codegeneration import expressions_template
from ffcx.codegeneration.backend import FFCXBackend
from ffcx.codegeneration.C.format_lines import format_indented_lines

logger = logging.getLogger("ffcx")


def generator(ir, parameters):
    """Generate UFC code for an expression."""

    logger.info("Generating code for expression:")
    logger.info("--- points: {}".format(ir.points))
    logger.info("--- name: {}".format(ir.name))

    factory_name = ir.name

    # Format declaration
    declaration = expressions_template.declaration.format(factory_name=factory_name)

    backend = FFCXBackend(ir, parameters)
    eg = ExpressionGenerator(ir, backend)

    code = {}
    code["name"] = "{}_expression".format(ir.name)
    parts = eg.generate()

    body = format_indented_lines(parts.cs_format(), 1)
    code["tabulate_expression"] = body

    code["original_coefficient_positions"] = format_indented_lines(
        eg.generate_original_coefficient_positions().cs_format(), 1)

    code["points"] = format_indented_lines(eg.generate_points().cs_format(), 1)

    if len(eg.ir.expression_shape) > 0:
        value_shape_decl = eg.generate_value_shape().cs_format()
    else:
        value_shape_decl = "static const int value_shape[1] = {0};"
    code["value_shape"] = format_indented_lines(value_shape_decl, 1)

    # Format implementation code
    implementation = expressions_template.factory.format(
        factory_name=factory_name,
        tabulate_expression=code["tabulate_expression"],
        original_coefficient_positions=code["original_coefficient_positions"],
        num_coefficients=len(ir.coefficient_numbering),
        num_points=ir.points.shape[0],
        topological_dimension=ir.points.shape[1],
        num_components=len(ir.expression_shape),
        points=code["points"],
        value_shape=code["value_shape"])

    return declaration, implementation


class ExpressionGenerator:
    def __init__(self, ir: dict, backend: FFCXBackend):

        if len(list(ir.integrand.keys())) != 1:
            raise RuntimeError("Only one set of points allowed for expression evaluation")

        self.ir = ir
        self.backend = backend
        self.scope = {}
        self._ufl_names = set()
        self.finalization_blocks = collections.defaultdict(list)
        self.symbol_counters = collections.defaultdict(int)
        self.shared_symbols = {}
        self.quadrature_rule = list(self.ir.integrand.keys())[0]

    def generate(self):
        L = self.backend.language

        parts = []

        parts += self.generate_element_tables()
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

    def generate_element_tables(self):
        """Generate tables of FE basis evaluated at specified points."""
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
        """Generate quadrature loop for this quadrature rule.

        In the context of expressions quadrature loop is not accumulated.

        """
        L = self.backend.language

        # Generate varying partition
        body = self.generate_varying_partition()
        body = L.commented_code_list(
            body, "Points loop body setup quadrature loop {}".format(self.quadrature_rule.id()))

        # Generate dofblock parts, some of this
        # will be placed before or after quadloop
        preparts, quadparts = \
            self.generate_dofblock_partition()
        body += quadparts

        # Wrap body in loop or scope
        if not body:
            # Could happen for integral with everything zero and optimized away
            quadparts = []
        else:
            iq = self.backend.symbols.quadrature_loop_index()
            num_points = self.quadrature_rule.points.shape[0]
            quadparts = [L.ForRange(iq, 0, num_points, body=body)]

        return preparts, quadparts

    def generate_varying_partition(self):
        """Generate factors of blocks which are not cellwise constant."""
        L = self.backend.language

        # Get annotated graph of factorisation
        F = self.ir.integrand[self.quadrature_rule]["factorization"]

        arraysymbol = L.Symbol("sv_{}".format(self.quadrature_rule.id()))
        parts = self.generate_partition(arraysymbol, F, "varying")
        parts = L.commented_code_list(
            parts, "Unstructured varying computations for quadrature rule {}".format(self.quadrature_rule.id()))
        return parts

    def generate_piecewise_partition(self):
        """Generate factors of blocks which are constant (i.e. do not depent on quadrature points)."""
        L = self.backend.language

        # Get annotated graph of factorisation
        F = self.ir.integrand[self.quadrature_rule]["factorization"]

        arraysymbol = L.Symbol("sp")
        parts = self.generate_partition(arraysymbol, F, "piecewise")
        parts = L.commented_code_list(parts, "Unstructured piecewise computations")
        return parts

    def generate_dofblock_partition(self):
        """Generate assignments of blocks multiplied with their factors into final tensor A."""
        block_contributions = self.ir.integrand[self.quadrature_rule]["block_contributions"]

        preparts = []
        quadparts = []

        blocks = [(blockmap, blockdata)
                  for blockmap, contributions in sorted(block_contributions.items())
                  for blockdata in contributions]

        for blockmap, blockdata in blocks:

            # Define code for block depending on mode
            block_preparts, block_quadparts = \
                self.generate_block_parts(blockmap, blockdata)

            # Add definitions
            preparts.extend(block_preparts)

            # Add computations
            quadparts.extend(block_quadparts)

        return preparts, quadparts

    def generate_block_parts(self, blockmap, blockdata):
        """Generate and return code parts for a given block."""
        L = self.backend.language

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
        A_shape = self.ir.tensor_shape
        Asym = self.backend.symbols.element_tensor()
        A = L.FlattenedArray(Asym, dims=[components] + [num_points] + A_shape)

        iq = self.backend.symbols.quadrature_loop_index()

        # Prepend dimensions of dofmap block with free index
        # for quadrature points and expression components
        B_indices = tuple([iq] + list(arg_indices))

        # Fetch code to access modified arguments
        # An access to FE table data
        arg_factors = self.get_arg_factors(blockdata, block_rank, B_indices)

        A_indices = []
        for i in range(len(blockmap)):
            offset = blockmap[i][0]
            A_indices.append(arg_indices[i] + offset)
        A_indices = tuple([iq] + A_indices)

        # Multiply collected factors
        # For each component of the factor expression
        # add result inside quadloop
        body = []

        for fi_ci in blockdata.factor_indices_comp_indices:
            f = self.get_var(F.nodes[fi_ci[0]]["expression"])
            Brhs = L.float_product([f] + arg_factors)
            body.append(L.AssignAdd(A[(fi_ci[1],) + A_indices], Brhs))

        for i in reversed(range(block_rank)):
            body = L.ForRange(
                B_indices[i + 1], 0, blockdims[i], body=body)
        quadparts += [body]

        return preparts, quadparts

    def get_arg_factors(self, blockdata, block_rank, indices):
        """Get argument factors (i.e. blocks).

        Parameters
        ----------
        blockdata
        block_rank
        indices
            Indices used to index element tables

        """
        L = self.backend.language

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
        L = self.backend.language
        name = "%s%d" % (basename, self.symbol_counters[basename])
        self.symbol_counters[basename] += 1
        return L.Symbol(name)

    def get_var(self, v):
        if v._ufl_is_literal_:
            return self.backend.ufl_to_language.get(v)
        f = self.scope.get(v)
        return f

    def generate_partition(self, symbol, F, mode):
        """Generate computations of factors of blocks."""
        L = self.backend.language

        definitions = []
        intermediates = []

        use_symbol_array = True

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
                vaccess = self.backend.access.get(mt.terminal, mt, tabledata, 0)
                vdef = self.backend.definitions.get(mt.terminal, mt, tabledata, 0, vaccess)

                # Store definitions of terminals in list
                assert isinstance(vdef, list)
                definitions.extend(vdef)
            else:
                # Get previously visited operands
                vops = [self.get_var(op) for op in v.ufl_operands]

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
                    if use_symbol_array:
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
            if use_symbol_array:
                alignas = self.ir.params["alignas"]
                parts += [L.ArrayDecl("ufc_scalar_t", symbol, len(intermediates), alignas=alignas)]
            parts += intermediates
        return parts

    def generate_original_coefficient_positions(self):
        """Generate original coefficient positions.

        Maps coefficient position index in processed expression
        to coefficient position index in original, non-processed expression.

        """
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
        """Generate and store compile-time known reference points at which the expression was evaluated."""
        L = self.backend.language
        parts = L.ArrayDecl("static const double", "points", values=self.ir.points,
                            sizes=self.ir.points.shape)
        return parts

    def generate_value_shape(self):
        L = self.backend.language
        parts = L.ArrayDecl("static const int", "value_shape", values=self.ir.expression_shape,
                            sizes=len(self.ir.expression_shape))
        return parts
