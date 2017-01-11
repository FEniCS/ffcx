# -*- coding: utf-8 -*-
# Copyright (C) 2011-2016 Martin Sandve Alnæs
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

"""Controlling algorithm for building the tabulate_tensor source structure from factorized representation."""

from collections import defaultdict

from ufl import product
from ufl.classes import ConstantValue, Condition
from ufl.measure import custom_integral_types, point_integral_types

from ffc.log import error, warning

from ffc.uflacs.build_uflacs_ir import get_common_block_data
from ffc.uflacs.elementtables import piecewise_ttypes

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

        # Cache of reusable blocks contributing to A
        self.shared_blocks = {}

        # Block contributions collected during generation to be added to A at the end
        self.finalization_blocks = defaultdict(list)

        # Set of counters used for assigning names to intermediate variables
        # TODO: Should this be part of the backend symbols? Doesn't really matter now.
        self.symbol_counters = defaultdict(int)


    def get_includes(self):
        "Return list of include statements needed to support generated code."
        includes = set()

        includes.add("#include <cstring>")  # for using memset
        #includes.add("#include <algorithm>")  # for using std::fill instead of memset

        cmath_names = set((
                "abs", "sign", "pow", "sqrt",
                "exp", "ln",
                "cos", "sin", "tan",
                "acos", "asin", "atan", "atan_2",
                "cosh", "sinh", "tanh",
                "acosh", "asinh", "atanh",
                "erf", "erfc",
            ))

        boost_math_names = set((
            "bessel_j", "bessel_y", "bessel_i", "bessel_k",
            ))

        # Only return the necessary headers
        if cmath_names & self._ufl_names:
            includes.add("#include <cmath>")

        if boost_math_names & self._ufl_names:
            includes.add("#include <boost/math/special_functions.hpp>")

        return sorted(includes)


    def init_scopes(self):
        "Initialize variable scope dicts."
        # Reset variables, separate sets for quadrature loop
        self.scopes = { num_points: {} for num_points in self.ir["all_num_points"] }
        self.scopes[None] = {}


    def set_var(self, num_points, v, vaccess):
        """"Set a new variable in variable scope dicts.

        Scope is determined by num_points which identifies the
        quadrature loop scope or None if outside quadrature loops.

        v is the ufl expression and vaccess is the CNodes
        expression to access the value in the code.
        """
        self.scopes[num_points][v] = vaccess


    def has_var(self, num_points, v):
        """"Check if variable exists in variable scope dicts.

        Return True if ufl expression v exists in the num_points scope.

        NB! Does not fall back to piecewise scope.
        """
        return v in self.scopes[num_points]


    def get_var(self, num_points, v):
        """"Lookup ufl expression v in variable scope dicts.

        Scope is determined by num_points which identifies the
        quadrature loop scope or None if outside quadrature loops.

        If v is not found in quadrature loop scope, the piecewise
        scope (None) is checked.

        Returns the CNodes expression to access the value in the code.
        """
        if v._ufl_is_literal_:
            return self.backend.ufl_to_language(v)
        f = self.scopes[num_points].get(v)
        if f is None:
            f = self.scopes[None][v]
        return f


    def new_temp_symbol(self, basename):
        "Create a new code symbol named basename + running counter."
        L = self.backend.language
        name = "%s%d" % (basename, self.symbol_counters[basename])
        self.symbol_counters[basename] += 1
        return L.Symbol(name)


    def get_temp_symbol(self, tempname, key):
        key = (tempname,) + key
        s = self.shared_blocks.get(key)
        defined = s is not None
        if not defined:
            s = self.new_temp_symbol(tempname)
            self.shared_blocks[key] = s
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

        # Generate code to set A = 0
        parts += self.generate_tensor_reset()

        # Generate code to compute piecewise constant scalar factors
        parts += self.generate_unstructured_piecewise_partition()

        # Loop generation code will produce parts to go before quadloops,
        # to define the quadloops, and to go after the quadloops
        all_preparts = []
        all_quadparts = []
        all_postparts = []

        # Go through each relevant quadrature loop
        if self.ir["integral_type"] in custom_integral_types:
            preparts, quadparts, postparts = \
              self.generate_runtime_quadrature_loop()
            all_preparts += preparts
            all_quadparts += quadparts
            all_postparts += postparts
        else:
            for num_points in self.ir["all_num_points"]:
                # Generate code to integrate reusable blocks of final element tensor
                preparts, quadparts, postparts = \
                    self.generate_quadrature_loop(num_points)
                all_preparts += preparts
                all_quadparts += quadparts
                all_postparts += postparts

        # Generate code to finish computing reusable blocks outside quadloop
        preparts, quadparts, postparts = \
            self.generate_dofblock_partition(None)
        all_preparts += preparts
        all_quadparts += quadparts
        all_postparts += postparts

        # Collect loop parts
        parts += all_preparts
        parts += all_quadparts
        parts += all_postparts

        # Generate code to add reusable blocks B* to element tensor A
        parts += self.generate_copyout_statements()

        return L.StatementList(parts)


    def generate_quadrature_tables(self):
        "Generate static tables of quadrature points and weights."
        L = self.backend.language

        parts = []

        # No quadrature tables for custom (given argument)
        # or point (evaluation in single vertex)
        skip = custom_integral_types + point_integral_types
        if self.ir["integral_type"] in skip:
            return parts

        # Loop over quadrature rules
        for num_points in self.ir["all_num_points"]:
            varying_ir = self.ir["varying_irs"][num_points]

            points, weights = self.ir["quadrature_rules"][num_points]
            assert num_points == len(weights)
            assert num_points == points.shape[0]

            # Generate quadrature weights array
            if varying_ir["need_weights"]:
                wsym = self.backend.symbols.weights_table(num_points)
                parts += [L.ArrayDecl("static const double", wsym, num_points, weights,
                                      alignas=self.ir["alignas"])]

            # Generate quadrature points array
            N = product(points.shape)
            if varying_ir["need_points"] and N:
                # Flatten array: (TODO: avoid flattening here, it makes padding harder)
                flattened_points = points.reshape(N)
                psym = self.backend.symbols.points_table(num_points)
                parts += [L.ArrayDecl("static const double", psym, N,
                                      flattened_points, alignas=self.ir["alignas"])]

        # Add leading comment if there are any tables
        parts = L.commented_code_list(parts, "Quadrature rules")
        return parts


    def generate_element_tables(self):
        """Generate static tables with precomputed element basis
        function values in quadrature points."""
        L = self.backend.language
        parts = []

        tables = self.ir["unique_tables"]
        table_types = self.ir["unique_table_types"]

        alignas = self.ir["alignas"]
        #padlen = self.ir["padlen"]

        if self.ir["integral_type"] in custom_integral_types:
            # Define only piecewise tables
            table_names = [name for name in sorted(tables)
                           if table_types[name] in piecewise_ttypes]
        else:
            # Define all tables
            table_names = sorted(tables)

        for name in table_names:
            table = tables[name]
            decl = L.ArrayDecl("static const double", name,
                               table.shape, table,
                               alignas=alignas)  # padlen=padlen)
            parts += [decl]

        # Add leading comment if there are any tables
        parts = L.commented_code_list(parts, [
            "Precomputed values of basis functions and precomputations",
            "FE* dimensions: [entities][points][dofs]",
            "PI* dimensions: [entities][dofs][dofs] or [entities][dofs]",
            "PM* dimensions: [entities][dofs][dofs]",
            ])
        return parts


    def generate_tensor_reset(self):
        "Generate statements for resetting the element tensor to zero."
        L = self.backend.language

        # TODO: Move this to language module, make CNode type
        def memzero(ptrname, size):
            tmp = "memset({ptrname}, 0, {size} * sizeof(*{ptrname}));"
            code = tmp.format(ptrname=str(ptrname), size=size)
            return L.VerbatimStatement(code)

        # Compute tensor size
        A = self.backend.symbols.element_tensor()
        A_size = product(self.ir["tensor_shape"])

        # Stitch it together
        parts = [L.Comment("Reset element tensor")]
        if A_size == 1:
            parts += [L.Assign(A[0], L.LiteralFloat(0.0))]
        else:
            parts += [memzero(A, A_size)]
        return parts


    def generate_quadrature_loop(self, num_points):
        "Generate quadrature loop with for this num_points."
        L = self.backend.language

        # Generate unstructured varying partition
        body = self.generate_unstructured_varying_partition(num_points)
        body = L.commented_code_list(body,
            "Quadrature loop body setup (num_points={0})".format(num_points))

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
        "Generate quadrature loop for custom integrals, with physical points given runtime."
        L = self.backend.language

        assert self.ir["integral_type"] in custom_integral_types

        num_points = self.ir["fake_num_points"]
        chunk_size = self.ir["chunk_size"]

        gdim = self.ir["geometric_dimension"]
        tdim = self.ir["topological_dimension"]

        alignas = self.ir["alignas"]
        # padlen = self.ir["padlen"]

        tables = self.ir["unique_tables"]
        table_types = self.ir["unique_table_types"]
        #table_origins = self.ir["unique_table_origins"]  # FIXME

        # Generate unstructured varying partition
        body = self.generate_unstructured_varying_partition(num_points)
        body = L.commented_code_list(body,
            ["Run-time quadrature loop body setup",
             "(chunk_size={0}, analysis_num_points={1})".format(
                 chunk_size, num_points)])

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


            ### Preparations for quadrature rules
            # 
            varying_ir = self.ir["varying_irs"][num_points]

            # Copy quadrature weights for this chunk
            if varying_ir["need_weights"]:
                cwsym = self.backend.symbols.custom_quadrature_weights()
                wsym = self.backend.symbols.custom_weights_table()
                rule_parts += [
                    L.ArrayDecl("double", wsym, chunk_size, 0,
                                alignas=alignas),
                    L.ForRange(iq, 0, num_points_in_block,
                               body=L.Assign(wsym[iq], cwsym[chunk_size*iq_chunk + iq])),
                    ]

            # Copy quadrature points for this chunk
            if varying_ir["need_points"]:
                cpsym = self.backend.symbols.custom_quadrature_points()
                psym = self.backend.symbols.custom_points_table()
                rule_parts += [
                    L.ArrayDecl("double", psym, chunk_size * gdim, 0,
                                alignas=alignas),
                    L.ForRange(iq, 0, num_points_in_block,
                        body=[L.Assign(psym[iq*gdim + i],
                                       cpsym[chunk_size*iq_chunk*gdim + iq*gdim + i])
                              for i in range(gdim)])
                    ]

            # Add leading comment if there are any tables
            rule_parts = L.commented_code_list(rule_parts,
                "Quadrature weights and points")


            ### Preparations for element tables
            table_parts = []

            # Only declare non-piecewise tables, computed inside chunk loop
            non_piecewise_tables = [name for name in sorted(tables)
                                    if table_types[name] not in piecewise_ttypes]
            for name in non_piecewise_tables:
                table = tables[name]
                decl = L.ArrayDecl("double", name,
                                   (1, chunk_size, table.shape[2]), 0,
                                   alignas=alignas)  # padlen=padlen)
                table_parts += [decl]

            table_parts += [L.Comment("FIXME: Fill element tables here")]
            #table_origins

            ### Gather all in chunk loop
            chunk_body = rule_parts + table_parts + [iq_body]
            quadparts = [L.ForRange(iq_chunk, 0, num_point_blocks, body=chunk_body)]

        return preparts, quadparts, postparts


    def generate_unstructured_piecewise_partition(self):
        L = self.backend.language

        num_points = None
        expr_ir = self.ir["piecewise_ir"]

        name = "sp"
        arraysymbol = L.Symbol(name)
        parts = self.generate_partition(arraysymbol,
                                        expr_ir["V"],
                                        expr_ir["V_active"],
                                        expr_ir["V_mts"],
                                        expr_ir["mt_tabledata"],
                                        num_points)
        parts = L.commented_code_list(parts,
            "Unstructured piecewise computations")
        return parts


    def generate_unstructured_varying_partition(self, num_points):
        L = self.backend.language

        expr_ir = self.ir["varying_irs"][num_points]

        name = "sv"
        arraysymbol = L.Symbol("%s%d" % (name, num_points))
        parts = self.generate_partition(arraysymbol,
                                        expr_ir["V"],
                                        expr_ir["V_varying"],
                                        expr_ir["V_mts"],
                                        expr_ir["mt_tabledata"],
                                        num_points)
        parts = L.commented_code_list(parts,
            "Unstructured varying computations for num_points=%d" % (num_points,))
        return parts


    def generate_partition(self, symbol, V, V_active, V_mts, mt_tabledata, num_points):
        L = self.backend.language

        definitions = []
        intermediates = []

        active_indices = [i for i, p in enumerate(V_active) if p]

        for i in active_indices:
            v = V[i]
            mt = V_mts[i]

            if v._ufl_is_literal_:
                vaccess = self.backend.ufl_to_language(v)
            elif mt is not None:
                # All finite element based terminals has table data, as well
                # as some but not all of the symbolic geometric terminals
                tabledata = mt_tabledata.get(mt)

                # Backend specific modified terminal translation
                vaccess = self.backend.access(mt.terminal,
                    mt, tabledata, num_points)
                vdef = self.backend.definitions(mt.terminal,
                    mt, tabledata, num_points, vaccess)

                # Store definitions of terminals in list
                assert isinstance(vdef, list)
                definitions.extend(vdef)
            else:
                # Get previously visited operands
                vops = [self.get_var(num_points, op) for op in v.ufl_operands]

                # Mapping UFL operator to target language
                self._ufl_names.add(v._ufl_handler_name_)
                vexpr = self.backend.ufl_to_language(v, *vops)

                # TODO: Let optimized ir provide mapping of vertex indices to
                # variable indices, marking which subexpressions to store in variables
                # and in what order:
                #j = variable_id[i]

                # Currently instead creating a new intermediate for
                # each subexpression except boolean conditions
                if isinstance(v, Condition):
                    # Inline the conditions x < y, condition values
                    # 'x' and 'y' may still be stored in intermediates.
                    # This removes the need to handle boolean intermediate variables.
                    # With tensor-valued conditionals it may not be optimal but we
                    # let the C++ compiler take responsibility for optimizing those cases.
                    j = None
                elif any(op._ufl_is_literal_ for op in v.ufl_operands):
                    # Skip intermediates for e.g. -2.0*x,
                    # resulting in lines like z = y + -2.0*x
                    j = None
                else:
                    j = len(intermediates)

                if j is not None:
                    # Record assignment of vexpr to intermediate variable
                    vaccess = symbol[j]
                    intermediates.append(L.Assign(vaccess, vexpr))
                else:
                    # Access the inlined expression
                    vaccess = vexpr

            # Store access node for future reference
            self.set_var(num_points, v, vaccess)

        # Join terminal computation, array of intermediate expressions,
        # and intermediate computations
        parts = []
        if definitions:
            parts += definitions
        if intermediates:
            parts += [L.ArrayDecl("double", symbol, len(intermediates),
                                  alignas=self.ir["alignas"])]
            parts += intermediates
        return parts


    def generate_dofblock_partition(self, num_points):
        L = self.backend.language

        if num_points is None:  # NB! None meaning piecewise partition, not custom integral
            block_contributions = self.ir["piecewise_ir"]["block_contributions"]
        else:
            block_contributions = self.ir["varying_irs"][num_points]["block_contributions"]

        preparts = []
        quadparts = []
        postparts = []

        for blockmap, contributions in sorted(block_contributions.items()):
            for blockdata in contributions:
                # Get symbol for already defined block B if it exists
                common_block_data = get_common_block_data(blockdata)
                B = self.shared_blocks.get(common_block_data)
                if B is None:
                    # Define code for block depending on mode
                    B, block_preparts, block_quadparts, block_postparts = \
                        self.generate_block_parts(num_points, blockmap, blockdata)

                    # Add definitions
                    preparts.extend(block_preparts)

                    # Add computations
                    quadparts.extend(block_quadparts)

                    # Add finalization
                    postparts.extend(block_postparts)

                    # Store reference for reuse
                    self.shared_blocks[common_block_data] = B

                # Add A[blockmap] += B[...] to finalization
                self.finalization_blocks[blockmap].append(B)

        return preparts, quadparts, postparts


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
        #tempnames = self.backend.symbols.block_temp_names
        #blocknames = self.backend.symbols.block_names
        tempnames = {
            #"preintegrated": "TI",
            "premultiplied": "TM",
            "partial": "TP",
            "full": "TF",
            "safe": "TS",
            "quadrature": "TQ",
            }
        blocknames = {
            #"preintegrated": "BI",
            #"premultiplied": "BM",
            #"partial": "BP",
            "full": "BF",
            "safe": "BS",
            "quadrature": "BQ",
            }

        fwtempname = "fw"
        tempname = tempnames.get(blockdata.block_mode)

        block_rank = len(blockmap)
        blockdims = tuple(len(dofmap) for dofmap in blockmap)

        ttypes = blockdata.ttypes
        if "zeros" in ttypes:
            error("Not expecting zero arguments to be left in dofblock generation.")

        if num_points is None:
            iq = None
        elif num_points == 1:
            iq = 0
        else:
            iq = self.backend.symbols.quadrature_loop_index()

        # Override dof index with quadrature loop index for arguments with
        # quadrature element, to index B like B[iq*num_dofs + iq]
        arg_indices = tuple(self.backend.symbols.argument_loop_index(i)
                            for i in range(block_rank))
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
            preparts.append(L.ArrayDecl("double", B, blockdims, 0,
                                        alignas=self.ir["alignas"]))

        # Get factor expression
        if blockdata.factor_is_piecewise:
            v = self.ir["piecewise_ir"]["V"][blockdata.factor_index]
        else:
            v = self.ir["varying_irs"][num_points]["V"][blockdata.factor_index]
        f = self.get_var(num_points, v)

        # Quadrature weight was removed in representation, add it back now
        if num_points is None:
            weight = L.LiteralFloat(1.0)
        elif self.ir["integral_type"] in custom_integral_types:
            weights = self.backend.symbols.custom_weights_table()
            weight = weights[iq]
        else:
            weights = self.backend.symbols.weights_table(num_points)
            weight = weights[iq]

        # Fetch code to access modified arguments
        if blockdata.block_mode in ("safe", "full", "partial"):
            arg_factors = []
            for i in range(block_rank):
                mad = blockdata.ma_data[i]
                td = mad.tabledata
                if td.is_piecewise:
                    scope = self.ir["piecewise_ir"]["modified_arguments"]
                else:
                    scope = self.ir["varying_irs"][num_points]["modified_arguments"]
                mt = scope[mad.ma_index]

                # Translate modified terminal to code
                # TODO: Move element table access out of backend?
                #       Not using self.backend.access.argument() here
                #       now because it assumes too much about indices.

                table = self.backend.symbols.element_table(td,
                    self.ir["entitytype"], mt.restriction)

                assert td.ttype != "zeros"

                if td.ttype == "ones":
                    arg_factor = L.LiteralFloat(1.0)
                elif td.ttype == "quadrature":  # TODO: Revisit all quadrature ttype checks
                    arg_factor = table[iq]
                else:
                    # Assuming B sparsity follows element table sparsity
                    arg_factor = table[B_indices[i]]
                arg_factors.append(arg_factor)
            assert len(arg_factors) == block_rank

        # Define fw = f * weight
        if blockdata.block_mode in ("safe", "full", "partial"):
            assert not blockdata.transposed, "Not handled yet"

            fw_rhs = L.float_product([f, weight])
            if not isinstance(fw_rhs, L.Product):
                fw = fw_rhs
            else:
                # Define and cache scalar temp variable
                key = (num_points, blockdata.factor_index, blockdata.factor_is_piecewise)
                fw, defined = self.get_temp_symbol(fwtempname, key)
                if not defined:
                    quadparts.append(L.VariableDecl("const double", fw, fw_rhs))

        if blockdata.block_mode == "safe":
            assert not blockdata.transposed

            # Collect scalar factors
            B_rhs = L.float_product([fw] + arg_factors)

            # Add result to block inside quadloop
            body = L.AssignAdd(B[B_indices], B_rhs)  # NB! += not =
            for i in reversed(range(block_rank)):
                if ttypes[i] != "quadrature":
                    body = L.ForRange(B_indices[i], 0, blockdims[i], body=body)
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
                    quadparts.append(L.ArrayDecl("double", P, P_dim, None,
                                                 alignas=self.ir["alignas"]))
                    P_rhs = L.float_product([fw, arg_factors[i]])
                    body = L.Assign(P[P_index], P_rhs)
                    #if ttypes[i] != "quadrature":  # FIXME: What does this mean here?
                    body = L.ForRange(P_index, 0, P_dim, body=body)
                    quadparts.append(body)

                B_rhs = P[P_index] * arg_factors[j]

            # Add result to block inside quadloop
            body = L.AssignAdd(B[B_indices], B_rhs)  # NB! += not =
            for i in reversed(range(block_rank)):
                if ttypes[i] != "quadrature":
                    body = L.ForRange(B_indices[i], 0, blockdims[i], body=body)

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
                preparts.append(L.ArrayDecl("double", P, P_dim, 0,
                                            alignas=self.ir["alignas"]))

                # Multiply collected factors
                P_rhs = L.float_product([fw, arg_factors[not_piecewise_index]])

                # Accumulate P += weight * f * args in quadrature loop
                body = L.AssignAdd(P[P_index], P_rhs)
                body = L.ForRange(P_index, 0, P_dim, body=body)
                quadparts.append(body)

            # Define B = B_rhs = piecewise_argument[:] * P[:], where P[:] = sum_q weight * f * other_argument[:]
            B_rhs = arg_factors[i] * P[P_index]

            # Define rhs expression for A[blockmap[arg_indices]] += A_rhs
            A_rhs = B_rhs

        elif blockdata.block_mode in ("premultiplied", "preintegrated"):
            if self.ir["integral_type"] == "interior_facet":
                # Get the facet entities
                entities = []
                for r in blockdata.restrictions:
                    if r is None:
                        entities.append(0)
                    else:
                        entities.append(self.backend.symbols.entity(self.ir["entitytype"], r))
                if blockdata.transposed:
                    P_block_indices = (entities[1], entities[0])
                else:
                    P_block_indices = tuple(entities)
            else:
                # Get the current cell or facet entity
                entity = self.backend.symbols.entity(self.ir["entitytype"], None)
                P_block_indices = (entity,)

            if blockdata.transposed:
                P_block_indices += (arg_indices[1], arg_indices[0])
            else:
                P_block_indices += arg_indices

            if blockdata.block_mode == "preintegrated":
                # Preintegrated should never get into quadloops
                assert num_points is None

                # Define B = B_rhs = f * PI where PI = sum_q weight * u * v
                PI = L.Symbol(blockdata.name)[P_block_indices]
                B_rhs = L.float_product([f, PI])

            elif blockdata.block_mode == "premultiplied":
                key = (num_points, blockdata.factor_index, blockdata.factor_is_piecewise)
                FI, defined = self.get_temp_symbol(tempname, key)
                if not defined:
                    # Declare FI = 0 before quadloop
                    preparts += [L.VariableDecl("double", FI, 0)]
                    # Accumulate FI += weight * f in quadparts
                    quadparts += [L.AssignAdd(FI, L.float_product([weight, f]))]

                # Define B_rhs = FI * PM where FI = sum_q weight*f, and PM = u * v
                PM = L.Symbol(blockdata.name)[P_block_indices]
                B_rhs = L.float_product([FI, PM])

            # Define rhs expression for A[blockmap[arg_indices]] += A_rhs
            A_rhs = B_rhs

        return A_rhs, preparts, quadparts, postparts


    def generate_copyout_statements(self):
        """Generate statements copying results to output array."""
        L = self.backend.language
        parts = []

        if self.ir["integral_type"] == "expression":
            # Not expecting any quadrature loop scopes here
            assert tuple(self.scopes.keys()) == (None,)

            # TODO: Get symbol from backend
            values = L.Symbol("values")

            # TODO: Allow expression compilation to compute multiple points at once!
            # Similarities to custom integrals in that points are given,
            # while different in output format: results are not accumulated
            # for each point but stored in output array instead.

            # Assign computed results to output variables
            pir = self.ir["piecewise_ir"]
            V = pir["V"]
            V_targets = pir["V_targets"]
            for i, fi in enumerate(V_targets):
                parts.append(L.Assign(values[i], self.get_var(None, V[fi])))
        else:
            # Get symbol, dimensions, and loop index symbols for A
            A = self.backend.symbols.element_tensor()
            A = L.FlattenedArray(A, dims=self.ir["tensor_shape"])
            rank = len(self.ir["tensor_shape"])
            indices = [self.backend.symbols.argument_loop_index(i)
                       for i in range(rank)]

            dofmaps = {}
            for blockmap, contributions in sorted(self.finalization_blocks.items()):
                # Define mapping from B indices to A indices
                A_indices = []
                for i in range(rank):
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
                            parts.append(L.ArrayDecl("static const int", DM, len(dofmap), dofmap))
                        j = DM[indices[i]]
                    A_indices.append(j)
                A_indices = tuple(A_indices)

                # TODO: If B[ij] = ... outside quadloop, change
                #       "A[ij] += B[ij]" into "A[ij] += ..." and drop the B[ij] temporary.

                # Sum up all blocks contributing to this blockmap
                if 0:
                    B_indices = indices
                    term = L.Sum([B[B_indices] for B in contributions])
                else: # XXX New suggestion, let B_indices be rolled into B_rhs
                    term = L.Sum([B_rhs for B_rhs in contributions])

                # TODO: need ttypes associated with this block to deal
                # with loop dropping for quadrature elements:
                ttypes = ()
                if ttypes == ("quadrature", "quadrature"):
                    debug("quadrature element block insertion not optimized")

                # Add components of all B's to A component in loop nest
                body = L.AssignAdd(A[A_indices], term)
                for i in reversed(range(rank)):
                    body = L.ForRange(indices[i], 0, len(blockmap[i]), body=body)

                # Add this block to parts
                parts.append(body)

        return parts
