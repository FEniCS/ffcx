# -*- coding: utf-8 -*-
# Copyright (C) 2013-2017 Martin Sandve Alnæs
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Controlling algorithm for building the tabulate_tensor
source structure from factorized representation."""

import collections
import itertools
import logging

import ufl
from ffc.codegeneration.backend import FFCBackend
from ffc.codegeneration.C.cnodes import pad_dim, pad_innermost_dim
from ffc.codegeneration.C.format_lines import format_indented_lines
from ffc.ir.representationutils import initialize_integral_code
from ffc.ir.uflacs.elementtables import piecewise_ttypes

logger = logging.getLogger(__name__)


def generate_integral_code(ir, parameters):
    """Generate code for integral from intermediate representation."""

    logger.info("Generating code from ffc.ir.uflacs representation")

    # FIXME: Is this the right precision value to use? Make it default to None or 0.
    precision = ir.integrals_metadata["precision"]

    # Create FFC C backend
    backend = FFCBackend(ir, parameters)

    # Configure kernel generator
    ig = IntegralGenerator(ir, backend, precision)

    # Generate code ast for the tabulate_tensor body
    parts = ig.generate()

    # Format code as string
    body = format_indented_lines(parts.cs_format(precision), 1)

    # Generate generic ffc code snippets and add uflacs specific parts
    code = initialize_integral_code(ir, parameters)
    code["tabulate_tensor"] = body

    return code


class IntegralGenerator(object):
    def __init__(self, ir, backend, precision):
        # Store ir
        self.ir = ir

        # Formatting precision
        self.precision = precision

        # Backend specific plugin with attributes
        # - language: for translating ufl operators to target language
        # - symbols: for translating ufl operators to target language
        # - definitions: for defining backend specific variables
        # - access: for accessing backend specific variables
        self.backend = backend

        # Set of operator names code has been generated for,
        # used in the end for selecting necessary includes
        self._ufl_names = set()

        # Initialize lookup tables for variable scopes
        self.init_scopes()

        # Cache
        self.shared_symbols = {}

        # Block contributions collected during generation to be added to A at the end
        self.finalization_blocks = collections.defaultdict(list)

        # Set of counters used for assigning names to intermediate variables
        self.symbol_counters = collections.defaultdict(int)

    def init_scopes(self):
        """Initialize variable scope dicts."""
        # Reset variables, separate sets for quadrature loop
        self.scopes = {num_points: {} for num_points in self.ir.all_num_points}
        self.scopes[None] = {}

    def set_var(self, num_points, v, vaccess):
        """Set a new variable in variable scope dicts.

        Scope is determined by num_points which identifies the
        quadrature loop scope or None if outside quadrature loops.

        v is the ufl expression and vaccess is the CNodes
        expression to access the value in the code.

        """
        self.scopes[num_points][v] = vaccess

    def get_var(self, num_points, v):
        """Lookup ufl expression v in variable scope dicts.

        Scope is determined by num_points which identifies the
        quadrature loop scope or None if outside quadrature loops.

        If v is not found in quadrature loop scope, the piecewise
        scope (None) is checked.

        Returns the CNodes expression to access the value in the code.
        """
        if v._ufl_is_literal_:
            return self.backend.ufl_to_language.get(v)
        f = self.scopes[num_points].get(v)
        if f is None:
            f = self.scopes[None][v]
        return f

    def new_temp_symbol(self, basename):
        """Create a new code symbol named basename + running counter."""
        L = self.backend.language
        name = "%s%d" % (basename, self.symbol_counters[basename])
        self.symbol_counters[basename] += 1
        return L.Symbol(name)

    def get_temp_symbol(self, tempname, key):
        key = (tempname, ) + key
        s = self.shared_symbols.get(key)
        defined = s is not None
        if not defined:
            s = self.new_temp_symbol(tempname)
            self.shared_symbols[key] = s
        return s, defined

    def generate(self):
        """Generate entire tabulate_tensor body.

        Assumes that the code returned from here will be wrapped in a context
        that matches a suitable version of the UFC tabulate_tensor signatures.
        """
        L = self.backend.language

        # Assert that scopes are empty: expecting this to be called only once
        assert not any(d for d in self.scopes.values())

        parts = []

        # Generate the tables of quadrature points and weights
        parts += self.generate_quadrature_tables()

        # Generate the tables of basis function values and preintegrated blocks
        parts += self.generate_element_tables()

        # Generate code to compute piecewise constant scalar factors
        parts += self.generate_unstructured_piecewise_partition()

        # Loop generation code will produce parts to go before quadloops,
        # to define the quadloops, and to go after the quadloops
        all_preparts = []
        all_quadparts = []
        all_postparts = []

        # Go through each relevant quadrature loop
        if self.ir.integral_type in ufl.measure.custom_integral_types:
            preparts, quadparts, postparts = \
                self.generate_runtime_quadrature_loop()
            all_preparts += preparts
            all_quadparts += quadparts
            all_postparts += postparts
        else:
            for num_points in self.ir.all_num_points:
                # Generate code to integrate reusable blocks of final element tensor
                preparts, quadparts, postparts = self.generate_quadrature_loop(num_points)
                all_preparts += preparts
                all_quadparts += quadparts
                all_postparts += postparts

        # Generate code to finish computing reusable blocks outside quadloop
        preparts, quadparts, postparts = \
            self.generate_dofblock_partition(None)
        all_preparts += preparts
        all_quadparts += quadparts
        all_postparts += postparts

        # Generate code to fill in A
        all_finalizeparts = []

        # Generate code to compute piecewise constant scalar factors
        # and set A at corresponding nonzero components
        all_finalizeparts += self.generate_preintegrated_dofblock_partition()

        # Generate code to add reusable blocks B* to element tensor A
        all_finalizeparts += self.generate_copyout_statements()

        # Collect parts before, during, and after quadrature loops
        parts += all_preparts
        parts += all_quadparts
        parts += all_postparts
        parts += all_finalizeparts

        return L.StatementList(parts)

    def generate_quadrature_tables(self):
        """Generate static tables of quadrature points and weights."""
        L = self.backend.language

        parts = []

        # No quadrature tables for custom (given argument)
        # or point (evaluation in single vertex)
        skip = ufl.measure.custom_integral_types + ufl.measure.point_integral_types
        if self.ir.integral_type in skip:
            return parts

        alignas = self.ir.params["alignas"]

        # Loop over quadrature rules
        for num_points in self.ir.all_num_points:
            varying_ir = self.ir.varying_irs[num_points]

            points, weights = self.ir.quadrature_rules[num_points]
            assert num_points == len(weights)
            assert num_points == points.shape[0]

            # Generate quadrature weights array
            if varying_ir["need_weights"]:
                wsym = self.backend.symbols.weights_table(num_points)
                parts += [
                    L.ArrayDecl(
                        "static const ufc_scalar_t", wsym, num_points, weights, alignas=alignas)
                ]

            # Generate quadrature points array
            N = ufl.product(points.shape)
            if varying_ir["need_points"] and N:
                # Flatten array: (TODO: avoid flattening here, it makes padding harder)
                flattened_points = points.reshape(N)
                psym = self.backend.symbols.points_table(num_points)
                parts += [
                    L.ArrayDecl(
                        "static const ufc_scalar_t", psym, N, flattened_points, alignas=alignas)
                ]

        # Add leading comment if there are any tables
        parts = L.commented_code_list(parts, "Quadrature rules")
        return parts

    def generate_element_tables(self):
        """Generate static tables with precomputed element basis
        function values in quadrature points."""
        L = self.backend.language
        parts = []

        tables = self.ir.unique_tables
        table_types = self.ir.unique_table_types
        inline_tables = self.ir.integral_type == "cell"

        alignas = self.ir.params["alignas"]
        padlen = self.ir.params["padlen"]

        if self.ir.integral_type in ufl.measure.custom_integral_types:
            # Define only piecewise tables
            table_names = [name for name in sorted(tables) if table_types[name] in piecewise_ttypes]
        else:
            # Define all tables
            table_names = sorted(tables)

        for name in table_names:
            table = tables[name]

            # Don't pad preintegrated tables
            if name[0] == "P":
                p = 1
            else:
                p = padlen

            # Skip tables that are inlined in code generation
            if inline_tables and name[:2] == "PI":
                continue

            decl = L.ArrayDecl(
                "static const ufc_scalar_t", name, table.shape, table, alignas=alignas, padlen=p)
            parts += [decl]

        # Add leading comment if there are any tables
        parts = L.commented_code_list(parts, [
            "Precomputed values of basis functions and precomputations",
            "FE* dimensions: [entities][points][dofs]",
            "PI* dimensions: [entities][dofs][dofs] or [entities][dofs]",
            "PM* dimensions: [entities][dofs][dofs]",
        ])
        return parts

    def generate_quadrature_loop(self, num_points):
        """Generate quadrature loop with for this num_points."""
        L = self.backend.language

        # Generate unstructured varying partition
        body = self.generate_unstructured_varying_partition(num_points)
        body = L.commented_code_list(
            body, "Quadrature loop body setup (num_points={0})".format(num_points))

        # Generate dofblock parts, some of this
        # will be placed before or after quadloop
        preparts, quadparts, postparts = \
            self.generate_dofblock_partition(num_points)
        body += quadparts

        # Wrap body in loop or scope
        if not body:
            # Could happen for integral with everything zero and optimized away
            quadparts = []
        elif num_points == 1:
            # For now wrapping body in Scope to avoid thinking about scoping issues
            quadparts = L.commented_code_list(L.Scope(body), "Only 1 quadrature point, no loop")
        else:
            # Regular case: define quadrature loop
            if num_points == 1:
                iq = 0
            else:
                iq = self.backend.symbols.quadrature_loop_index()
            quadparts = [L.ForRange(iq, 0, num_points, body=body)]

        return preparts, quadparts, postparts

    def generate_runtime_quadrature_loop(self):
        """Generate quadrature loop for custom integrals, with physical points given runtime."""
        L = self.backend.language

        assert self.ir.integral_type in ufl.measure.custom_integral_types

        num_points = self.ir.fake_num_points
        chunk_size = self.ir.params["chunk_size"]

        gdim = self.ir.geometric_dimension

        alignas = self.ir.params["alignas"]

        tables = self.ir.unique_tables
        table_types = self.ir.unique_table_types

        # Generate unstructured varying partition
        body = self.generate_unstructured_varying_partition(num_points)
        body = L.commented_code_list(body, [
            "Run-time quadrature loop body setup",
            "(chunk_size={0}, analysis_num_points={1})".format(chunk_size, num_points)
        ])

        # Generate dofblock parts, some of this
        # will be placed before or after quadloop
        preparts, quadparts, postparts = \
            self.generate_dofblock_partition(num_points)
        body += quadparts

        # Wrap body in loop
        if not body:
            # Could happen for integral with everything zero and optimized away
            quadparts = []
        else:
            rule_parts = []

            # Define two-level quadrature loop; over chunks then over points in chunk
            iq_chunk = L.Symbol("iq_chunk")
            np = self.backend.symbols.num_custom_quadrature_points()
            num_point_blocks = (np + chunk_size - 1) / chunk_size
            iq = self.backend.symbols.quadrature_loop_index()

            # Not assuming runtime size to be multiple by chunk size
            num_points_in_block = L.Symbol("num_points_in_chunk")
            decl = L.VariableDecl("const int", num_points_in_block,
                                  L.Call("min", (chunk_size, np - iq_chunk * chunk_size)))
            rule_parts.append(decl)

            iq_body = L.ForRange(iq, 0, num_points_in_block, body=body)

            # Preparations for quadrature rules
            #
            varying_ir = self.ir.varying_irs[num_points]

            # Copy quadrature weights for this chunk
            if varying_ir["need_weights"]:
                cwsym = self.backend.symbols.custom_quadrature_weights()
                wsym = self.backend.symbols.custom_weights_table()
                rule_parts += [
                    L.ArrayDecl("ufc_scalar_t", wsym, chunk_size, 0, alignas=alignas),
                    L.ForRange(
                        iq,
                        0,
                        num_points_in_block,
                        body=L.Assign(wsym[iq], cwsym[chunk_size * iq_chunk + iq])),
                ]

            # Copy quadrature points for this chunk
            if varying_ir["need_points"]:
                cpsym = self.backend.symbols.custom_quadrature_points()
                psym = self.backend.symbols.custom_points_table()
                rule_parts += [
                    L.ArrayDecl("ufc_scalar_t", psym, chunk_size * gdim, 0, alignas=alignas),
                    L.ForRange(
                        iq,
                        0,
                        num_points_in_block,
                        body=[
                            L.Assign(psym[iq * gdim + i],
                                     cpsym[chunk_size * iq_chunk * gdim + iq * gdim + i])
                            for i in range(gdim)
                        ])
                ]

            # Add leading comment if there are any tables
            rule_parts = L.commented_code_list(rule_parts, "Quadrature weights and points")

            # Preparations for element tables
            table_parts = []

            # Only declare non-piecewise tables, computed inside chunk loop
            non_piecewise_tables = [
                name for name in sorted(tables) if table_types[name] not in piecewise_ttypes
            ]
            for name in non_piecewise_tables:
                table = tables[name]
                decl = L.ArrayDecl(
                    "ufc_scalar_t", name, (1, chunk_size, table.shape[2]), 0,
                    alignas=alignas)  # padlen=padlen)
                table_parts += [decl]

            table_parts += [L.Comment("FIXME: Fill element tables here")]
            # table_origins

            # Gather all in chunk loop
            chunk_body = rule_parts + table_parts + [iq_body]
            quadparts = [L.ForRange(iq_chunk, 0, num_point_blocks, body=chunk_body)]

        return preparts, quadparts, postparts

    def generate_unstructured_piecewise_partition(self):
        L = self.backend.language

        # Get annotated graph of factorisation
        F = self.ir.piecewise_ir["factorization"]

        arraysymbol = L.Symbol("sp")
        num_points = None
        parts = self.generate_partition(arraysymbol, F, "piecewise", num_points)
        parts = L.commented_code_list(parts, "Unstructured piecewise computations")
        return parts

    def generate_unstructured_varying_partition(self, num_points):
        L = self.backend.language

        # Get annotated graph of factorisation
        F = self.ir.varying_irs[num_points]["factorization"]

        arraysymbol = L.Symbol("sv%d" % num_points)
        parts = self.generate_partition(arraysymbol, F, "varying", num_points)
        parts = L.commented_code_list(parts, "Unstructured varying computations for num_points=%d" %
                                      (num_points, ))
        return parts

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
            self.set_var(num_points, v, vaccess)

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

    def generate_dofblock_partition(self, num_points):
        if num_points is None:  # NB! None meaning piecewise partition, not custom integral
            block_contributions = self.ir.piecewise_ir["block_contributions"]
        else:
            block_contributions = self.ir.varying_irs[num_points]["block_contributions"]

        preparts = []
        quadparts = []
        postparts = []

        blocks = [(blockmap, blockdata)
                  for blockmap, contributions in sorted(block_contributions.items())
                  for blockdata in contributions if blockdata.block_mode != "preintegrated"]

        for blockmap, blockdata in blocks:

            # Define code for block depending on mode
            B, block_preparts, block_quadparts, block_postparts = \
                self.generate_block_parts(num_points, blockmap, blockdata)

            # Add definitions
            preparts.extend(block_preparts)

            # Add computations
            quadparts.extend(block_quadparts)

            # Add finalization
            postparts.extend(block_postparts)

            # Add A[blockmap] += B[...] to finalization
            self.finalization_blocks[blockmap].append(B)

        return preparts, quadparts, postparts

    def get_entities(self, blockdata):
        L = self.backend.language

        if self.ir.integral_type == "interior_facet":
            # Get the facet entities
            entities = []
            for r in blockdata.restrictions:
                if r is None:
                    entities.append(0)
                else:
                    entities.append(self.backend.symbols.entity(self.ir.entitytype, r))
            if blockdata.transposed:
                return (entities[1], entities[0])
            else:
                return tuple(entities)
        else:
            # Get the current cell or facet entity
            if blockdata.is_uniform:
                # uniform, i.e. constant across facets
                entity = L.LiteralInt(0)
            else:
                entity = self.backend.symbols.entity(self.ir.entitytype, None)
            return (entity, )

    def get_arg_factors(self, blockdata, block_rank, num_points, iq, indices):
        L = self.backend.language

        arg_factors = []
        for i in range(block_rank):
            mad = blockdata.ma_data[i]
            td = mad.tabledata
            if td.is_piecewise:
                scope = self.ir.piecewise_ir["modified_arguments"]
            else:
                scope = self.ir.varying_irs[num_points]["modified_arguments"]
            mt = scope[mad.ma_index]

            # Translate modified terminal to code
            # TODO: Move element table access out of backend?
            #       Not using self.backend.access.argument() here
            #       now because it assumes too much about indices.

            table = self.backend.symbols.element_table(td, self.ir.entitytype, mt.restriction)

            assert td.ttype != "zeros"

            if td.ttype == "ones":
                arg_factor = L.LiteralFloat(1.0)
            elif td.ttype == "quadrature":  # TODO: Revisit all quadrature ttype checks
                arg_factor = table[iq]
            else:
                # Assuming B sparsity follows element table sparsity
                arg_factor = table[indices[i]]
            arg_factors.append(arg_factor)
        return arg_factors

    def generate_block_parts(self, num_points, blockmap, blockdata):
        """Generate and return code parts for a given block.

        Returns parts occuring before, inside, and after
        the quadrature loop identified by num_points.

        Should be called with num_points=None for quadloop-independent blocks.
        """
        L = self.backend.language

        # The parts to return
        preparts = []
        quadparts = []
        postparts = []

        # TODO: Define names in backend symbols?
        # tempnames = self.backend.symbols.block_temp_names
        # blocknames = self.backend.symbols.block_names
        tempnames = {
            # "preintegrated": "TI",
            "premultiplied": "TM",
            "partial": "TP",
            "full": "TF",
            "safe": "TS",
            "quadrature": "TQ",
        }
        blocknames = {
            # "preintegrated": "BI",
            # "premultiplied": "BM",
            # "partial": "BP",
            "full": "BF",
            "safe": "BS",
            "quadrature": "BQ",
        }

        tempname = tempnames.get(blockdata.block_mode)

        alignas = self.ir.params["alignas"]
        padlen = self.ir.params["padlen"]

        block_rank = len(blockmap)
        blockdims = tuple(len(dofmap) for dofmap in blockmap)
        padded_blockdims = pad_innermost_dim(blockdims, padlen)

        ttypes = blockdata.ttypes
        if "zeros" in ttypes:
            raise RuntimeError("Not expecting zero arguments to be left in dofblock generation.")

        if num_points is None:
            iq = None
        elif num_points == 1:
            iq = 0
        else:
            iq = self.backend.symbols.quadrature_loop_index()

        # Override dof index with quadrature loop index for arguments with
        # quadrature element, to index B like B[iq*num_dofs + iq]
        arg_indices = tuple(self.backend.symbols.argument_loop_index(i) for i in range(block_rank))
        B_indices = []
        for i in range(block_rank):
            if ttypes[i] == "quadrature":
                B_indices.append(iq)
            else:
                B_indices.append(arg_indices[i])
        B_indices = tuple(B_indices)

        # Define unique block symbol
        blockname = blocknames.get(blockdata.block_mode)
        if blockname:
            B = self.new_temp_symbol(blockname)
            # Add initialization of this block to parts
            # For all modes, block definition occurs before quadloop
            preparts.append(
                L.ArrayDecl("ufc_scalar_t", B, blockdims, 0, alignas=alignas, padlen=padlen))

        # Get factor expression
        if blockdata.factor_is_piecewise:
            F = self.ir.piecewise_ir["factorization"]
        else:
            F = self.ir.varying_irs[num_points]["factorization"]
        v = F.nodes[blockdata.factor_index]['expression']
        f = self.get_var(num_points, v)

        # Quadrature weight was removed in representation, add it back now
        if num_points is None:
            weight = L.LiteralFloat(1.0)
        elif self.ir.integral_type in ufl.measure.custom_integral_types:
            weights = self.backend.symbols.custom_weights_table()
            weight = weights[iq]
        else:
            weights = self.backend.symbols.weights_table(num_points)
            weight = weights[iq]

        # Define fw = f * weight
        if blockdata.block_mode in ("safe", "full", "partial"):
            assert not blockdata.transposed, "Not handled yet"

            # Fetch code to access modified arguments
            arg_factors = self.get_arg_factors(blockdata, block_rank, num_points, iq, B_indices)

            fw_rhs = L.float_product([f, weight])
            if not isinstance(fw_rhs, L.Product):
                fw = fw_rhs
            else:
                # Define and cache scalar temp variable
                key = (num_points, blockdata.factor_index, blockdata.factor_is_piecewise)
                fw, defined = self.get_temp_symbol("fw", key)
                if not defined:
                    quadparts.append(L.VariableDecl("const ufc_scalar_t", fw, fw_rhs))

                # Plan for vectorization of fw computations over iq:
                # 1) Define fw as arrays e.g. "double fw0[nq];" outside quadloop
                # 2) Access as fw0[iq] of course
                # 3) Split quadrature loops, one for fw computation and one for blocks
                # 4) Pad quadrature rule with 0 weights and last point

                # Plan for vectorization of coefficient evaluation over iq:
                # 1) Define w0_c1 etc as arrays e.g. "double w0_c1[nq] = {};" outside quadloop
                # 2) Access as w0_c1[iq] of course
                # 3) Splitquadrature loops, coefficients before fw computation
                # 4) Possibly swap loops over iq and ic:
                #    for(ic) for(iq) w0_c1[iq] = w[0][ic] * FE[iq][ic];

        if blockdata.block_mode == "safe":
            # Naively accumulate integrand for this block in the innermost loop
            assert not blockdata.transposed
            B_rhs = L.float_product([fw] + arg_factors)
            body = L.AssignAdd(B[B_indices], B_rhs)  # NB! += not =
            for i in reversed(range(block_rank)):
                body = L.ForRange(B_indices[i], 0, padded_blockdims[i], body=body)
            quadparts += [body]

            # Define rhs expression for A[blockmap[arg_indices]] += A_rhs
            A_rhs = B[arg_indices]

        elif blockdata.block_mode == "full":
            assert not blockdata.transposed, "Not handled yet"

            if block_rank < 2:
                # Multiply collected factors
                B_rhs = L.float_product([fw] + arg_factors)
            else:
                # TODO: Pick arg with smallest dimension, or pick
                # based on global optimization to reuse more blocks
                i = 0  # Index selected for precomputation
                j = 1 - i

                P_index = B_indices[i]

                key = (num_points, blockdata.factor_index, blockdata.factor_is_piecewise,
                       arg_factors[i].ce_format(self.precision))
                P, defined = self.get_temp_symbol(tempname, key)
                if not defined:
                    # TODO: If FE table is varying and only used in contexts
                    # where it's multiplied by weight, we can premultiply it!
                    # Then this would become P = f * preweighted_FE_table[:].

                    # Define and compute intermediate value
                    # P[:] = (weight * f) * args[i][:]
                    # inside quadrature loop
                    P_dim = blockdims[i]
                    quadparts.append(
                        L.ArrayDecl("ufc_scalar_t", P, P_dim, None, alignas=alignas, padlen=padlen))
                    P_rhs = L.float_product([fw, arg_factors[i]])
                    body = L.Assign(P[P_index], P_rhs)
                    # if ttypes[i] != "quadrature":  # FIXME: What does this mean here?
                    vectorize = self.ir.params["vectorize"]
                    body = L.ForRange(P_index, 0, P_dim, body=body, vectorize=vectorize)
                    quadparts.append(body)

                B_rhs = P[P_index] * arg_factors[j]

            # Add result to block inside quadloop
            body = L.AssignAdd(B[B_indices], B_rhs)  # NB! += not =
            for i in reversed(range(block_rank)):
                # Vectorize only the innermost loop
                vectorize = self.ir.params["vectorize"] and (i == block_rank - 1)
                if ttypes[i] != "quadrature":
                    body = L.ForRange(
                        B_indices[i], 0, padded_blockdims[i], body=body, vectorize=vectorize)
            quadparts += [body]

            # Define rhs expression for A[blockmap[arg_indices]] += A_rhs
            A_rhs = B[arg_indices]

        elif blockdata.block_mode == "partial":
            # TODO: To handle transpose here, must add back intermediate block B
            assert not blockdata.transposed, "Not handled yet"

            # Get indices and dimensions right here...
            assert block_rank == 2
            i = blockdata.piecewise_ma_index
            not_piecewise_index = 1 - i

            P_index = arg_indices[not_piecewise_index]

            key = (num_points, blockdata.factor_index, blockdata.factor_is_piecewise,
                   arg_factors[not_piecewise_index].ce_format(self.precision))
            P, defined = self.get_temp_symbol(tempname, key)
            if not defined:
                # Declare P table in preparts
                P_dim = blockdims[not_piecewise_index]
                preparts.append(
                    L.ArrayDecl("ufc_scalar_t", P, P_dim, 0, alignas=alignas, padlen=padlen))

                # Multiply collected factors
                P_rhs = L.float_product([fw, arg_factors[not_piecewise_index]])

                # Accumulate P += weight * f * args in quadrature loop
                body = L.AssignAdd(P[P_index], P_rhs)
                body = L.ForRange(P_index, 0, pad_dim(P_dim, padlen), body=body)
                quadparts.append(body)

            # Define B = B_rhs = piecewise_argument[:] * P[:],
            # where P[:] = sum_q weight * f * other_argument[:]
            B_rhs = arg_factors[i] * P[P_index]

            # Define rhs expression for A[blockmap[arg_indices]] += A_rhs
            A_rhs = B_rhs

        elif blockdata.block_mode in ("premultiplied", "preintegrated"):
            P_ii = self.get_entities(blockdata)
            if blockdata.transposed:
                P_ii += arg_indices[::-1]
            else:
                P_ii += arg_indices

            if blockdata.block_mode == "preintegrated":
                # Preintegrated should never get into quadloops
                assert num_points is None

                # Define B = B_rhs = f * PI where PI = sum_q weight * u * v
                PI = L.Symbol(blockdata.name)[P_ii]
                B_rhs = L.float_product([f, PI])

            elif blockdata.block_mode == "premultiplied":
                key = (num_points, blockdata.factor_index, blockdata.factor_is_piecewise)
                FI, defined = self.get_temp_symbol(tempname, key)
                if not defined:
                    # Declare FI = 0 before quadloop
                    preparts += [L.VariableDecl("ufc_scalar_t", FI, 0)]
                    # Accumulate FI += weight * f in quadparts
                    quadparts += [L.AssignAdd(FI, L.float_product([weight, f]))]

                # Define B_rhs = FI * PM where FI = sum_q weight*f, and PM = u * v
                PM = L.Symbol(blockdata.name)[P_ii]
                B_rhs = L.float_product([FI, PM])

            # Define rhs expression for A[blockmap[arg_indices]] += A_rhs
            A_rhs = B_rhs

        # Equip code with comments
        comments = ["UFLACS block mode: {}".format(blockdata.block_mode)]
        preparts = L.commented_code_list(preparts, comments)
        quadparts = L.commented_code_list(quadparts, comments)
        postparts = L.commented_code_list(postparts, comments)

        return A_rhs, preparts, quadparts, postparts

    def generate_preintegrated_dofblock_partition(self):
        # FIXME: Generalize this to unrolling all A[] += ... loops,
        # or all loops with noncontiguous DM??
        L = self.backend.language

        block_contributions = self.ir.piecewise_ir["block_contributions"]

        blocks = [(blockmap, blockdata)
                  for blockmap, contributions in sorted(block_contributions.items())
                  for blockdata in contributions if blockdata.block_mode == "preintegrated"]

        # Get symbol, dimensions, and loop index symbols for A
        A_shape = self.ir.tensor_shape
        A_size = ufl.product(A_shape)
        A_rank = len(A_shape)

        # TODO: there's something like shape2strides(A_shape) somewhere
        # A_strides = ufl.utils.indexflattening.shape_to_strides(A_shape)

        A_strides = [1] * A_rank
        for i in reversed(range(0, A_rank - 1)):
            A_strides[i] = A_strides[i + 1] * A_shape[i + 1]

        A_values = [0.0] * A_size

        for blockmap, blockdata in blocks:
            # Accumulate A[blockmap[...]] += f*PI[...]

            # Get table for inlining
            tables = self.ir.unique_tables
            table = tables[blockdata.name]
            inline_table = self.ir.integral_type == "cell"

            # Get factor expression
            v = self.ir.piecewise_ir["factorization"].nodes[blockdata.factor_index]['expression']
            f = self.get_var(None, v)

            # Define rhs expression for A[blockmap[arg_indices]] += A_rhs
            # A_rhs = f * PI where PI = sum_q weight * u * v
            PI = L.Symbol(blockdata.name)

            # Define indices into preintegrated block
            P_entity_indices = self.get_entities(blockdata)
            if inline_table:
                assert P_entity_indices == (L.LiteralInt(0), )
                assert table.shape[0] == 1

            # Unroll loop
            blockshape = [len(DM) for DM in blockmap]
            blockrange = [range(d) for d in blockshape]

            for ii in itertools.product(*blockrange):
                A_ii = sum(A_strides[i] * blockmap[i][ii[i]] for i in range(len(ii)))
                if blockdata.transposed:
                    P_arg_indices = (ii[1], ii[0])
                else:
                    P_arg_indices = ii

                if inline_table:
                    # Extract float value of PI[P_ii]
                    Pval = table[0]  # always entity 0
                    for i in P_arg_indices:
                        Pval = Pval[i]
                    A_rhs = Pval * f
                else:
                    # Index the static preintegrated table:
                    P_ii = P_entity_indices + P_arg_indices
                    A_rhs = f * PI[P_ii]

                A_values[A_ii] = A_values[A_ii] + A_rhs

        code = self.generate_tensor_value_initialization(A_values)
        return L.commented_code_list(code, "UFLACS block mode: preintegrated")

    def generate_tensor_value_initialization(self, A_values):
        parts = []

        L = self.backend.language
        A = self.backend.symbols.element_tensor()
        A_size = len(A_values)

        init_mode = self.ir.params["tensor_init_mode"]
        z = L.LiteralFloat(0.0)

        k = L.Symbol("k")  # Index for zeroing arrays

        if init_mode == "direct":
            # Generate A[i] = A_values[i] including zeros
            for i in range(A_size):
                parts += [L.Assign(A[i], A_values[i])]
        elif init_mode == "upfront":
            # Zero everything first
            parts += [L.ForRange(k, 0, A_size, index_type="int", body=L.Assign(A[k], 0.0))]

            # Generate A[i] = A_values[i] skipping zeros
            for i in range(A_size):
                if not (A_values[i] == 0.0 or A_values[i] == z):
                    parts += [L.Assign(A[i], A_values[i])]
        elif init_mode == "interleaved":
            # Generate A[i] = A_values[i] with interleaved zero filling
            i = 0
            zero_begin = 0
            zero_end = zero_begin
            while i < A_size:
                if A_values[i] == 0.0 or A_values[i] == z:
                    # Update range of A zeros
                    zero_end = i + 1
                else:
                    # Set zeros of A just prior to A[i]
                    if zero_end == zero_begin + 1:
                        parts += [L.Assign(A[zero_begin], 0.0)]
                    elif zero_end > zero_begin:
                        parts += [
                            L.ForRange(
                                k, zero_begin, zero_end, index_type="int", body=L.Assign(A[k], 0.0))
                        ]
                    zero_begin = i + 1
                    zero_end = zero_begin
                    # Set A[i] value
                    parts += [L.Assign(A[i], A_values[i])]
                i += 1
            if zero_end == zero_begin + 1:
                parts += [L.Assign(A[zero_begin], 0.0)]
            elif zero_end > zero_begin:
                parts += [
                    L.ForRange(k, zero_begin, zero_end, index_type="int", body=L.Assign(A[k], 0.0))
                ]
        else:
            raise RuntimeError("Invalid init_mode parameter %s" % (init_mode, ))

        return parts

    def generate_copyout_statements(self):
        L = self.backend.language
        parts = []

        # Get symbol, dimensions, and loop index symbols for A
        A_shape = self.ir.tensor_shape
        A_rank = len(A_shape)

        Asym = self.backend.symbols.element_tensor()
        A = L.FlattenedArray(Asym, dims=A_shape)

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
            A_indices = tuple(A_indices)

            # Sum up all blocks contributing to this blockmap
            term = L.Sum([B_rhs for B_rhs in contributions])

            # TODO: need ttypes associated with this block to deal
            # with loop dropping for quadrature elements:
            ttypes = ()
            if ttypes == ("quadrature", "quadrature"):
                logger.debug("quadrature element block insertion not optimized")

            # Add components of all B's to A component in loop nest
            body = L.AssignAdd(A[A_indices], term)
            for i in reversed(range(A_rank)):
                body = L.ForRange(indices[i], 0, len(blockmap[i]), body=body)

            # Add this block to parts
            parts.append(body)

        # Place static dofmap tables first
        parts = dofmap_parts + parts

        return parts
