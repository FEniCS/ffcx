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

from ffc.uflacs.analysis.modified_terminals import analyse_modified_terminal, is_modified_terminal


class IntegralGenerator(object):

    def __init__(self, ir, backend):
        # Store ir
        self.ir = ir

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


    def generate(self):
        """Generate entire tabulate_tensor body.

        Assumes that the code returned from here will be wrapped in a context
        that matches a suitable version of the UFC tabulate_tensor signatures.
        """
        L = self.backend.language

        # Assert that scopes are empty: expecting this to be called only once
        assert not any(d for d in self.scopes.values())

        parts = []
        parts += self.generate_quadrature_tables()
        parts += self.generate_element_tables()
        parts += self.generate_tensor_reset()
        parts += self.generate_unstructured_piecewise_partition()

        all_preparts = []
        all_quadparts = []
        all_postparts = []

        for num_points in self.ir["all_num_points"]:
            preparts, quadparts, postparts = \
                self.generate_quadrature_loop(num_points)
            all_preparts += preparts
            all_quadparts += quadparts
            all_postparts += postparts

        preparts, quadparts, postparts = self.generate_dofblock_partition(None)
        parts += all_preparts
        parts += preparts
        parts += all_quadparts
        parts += quadparts
        parts += all_postparts
        parts += postparts

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
                wsym = self.backend.symbols.weights_array(num_points)
                parts += [L.ArrayDecl("static const double", wsym, num_points, weights,
                                      alignas=self.ir["alignas"])]

            # Generate quadrature points array
            N = product(points.shape)
            if varying_ir["need_points"] and N:
                # Flatten array: (TODO: avoid flattening here, it makes padding harder)
                flattened_points = points.reshape(N)
                psym = self.backend.symbols.points_array(num_points)
                parts += [L.ArrayDecl("static const double", psym, N,
                                      flattened_points, alignas=self.ir["alignas"])]

        # Add leading comment if there are any tables
        parts = L.commented_code_list(parts,
            "Section for quadrature weights and points")
        return parts


    def generate_element_tables(self):
        """Generate static tables with precomputed element basis
        function values in quadrature points."""
        L = self.backend.language
        parts = []

        tables = self.ir["unique_tables"]
        for name in sorted(tables):
            # TODO: Not padding, consider when and if to do so
            table = tables[name]
            decl = L.ArrayDecl("static const double", name,
                               table.shape, table,
                               alignas=self.ir["alignas"])
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
            iq = self.backend.symbols.quadrature_loop_index(num_points)
            np = self.backend.symbols.num_quadrature_points(num_points)
            quadparts = [L.ForRange(iq, 0, np, body=body)]

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
                                        expr_ir["V_table_data"],
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
                                        expr_ir["V_table_data"],
                                        num_points)
        parts = L.commented_code_list(parts,
            "Unstructured varying computations for num_points=%d" % (num_points,))
        return parts


    def generate_partition(self, symbol, V, V_active, V_table_data, num_points):
        L = self.backend.language

        definitions = []
        intermediates = []

        active_indices = [i for i, p in enumerate(V_active) if p]

        for i in active_indices:
            v = V[i]

            # XXX: Enable this after tests pass again to avoid too much changes at once:
            #if v._ufl_is_literal_:
            #    vaccess = self.backend.ufl_to_language(v)
            #elif is_modified_terminal(v):
            if is_modified_terminal(v):
                mt = analyse_modified_terminal(v)

                tabledata = V_table_data[i]
                #tabledata = modified_terminal_table_data[mt]

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
                # XXX: Enable this after tests pass again to avoid too much changes at once:
                #elif any(op._ufl_is_literal_ for op in v.ufl_operands):
                #    # Skip intermediates for e.g. -2.0*x,
                #    # resulting in lines like z = y + -2.0*x
                #    j = None
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

        for dofblock, contributions in sorted(block_contributions.items()):
            blockdims = tuple(end - begin for begin, end in dofblock)
            for blockdata in contributions:
                # Get symbol for already defined block B if it exists
                B = self.shared_blocks.get(blockdata)
                if B is None:
                    # Define code for block depending on mode
                    B, block_preparts, block_quadparts, block_postparts = \
                        self.generate_block_parts(num_points, dofblock, blockdata)

                    # Add definitions
                    preparts.extend(block_preparts)

                    # Add computations
                    quadparts.extend(block_quadparts)

                    # Add finalization
                    postparts.extend(block_postparts)

                    # Store reference for reuse
                    self.shared_blocks[blockdata] = B

                # Add A[dofblock] += B[...] to finalization
                self.finalization_blocks[dofblock].append(B)

        return preparts, quadparts, postparts


    def generate_block_parts(self, num_points, dofblock, blockdata):
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
        blocknames = {
            "preintegrated": "BI",
            "premultiplied": "BM",
            "partial": "BP",
            "full": "BF",
            }

        block_rank = len(dofblock)
        blockdims = tuple(end - begin for begin, end in dofblock)
        ttypes = blockdata.ttypes
        if "zeros" in ttypes:
            error("Not expecting zero arguments to be left in dofblock generation.")

        if num_points is not None:
            iq = self.backend.symbols.quadrature_loop_index(num_points)

        # Override dof index with quadrature loop index for arguments with
        # quadrature element, to index B like B[iq*num_dofs + iq]
        B_indices = []
        for i in range(block_rank):
            if ttypes[i] == "quadrature":
                B_indices.append(iq)
            else:
                B_indices.append(self.backend.symbols.argument_loop_index(i))
        B_indices = tuple(B_indices)

        # Define unique block symbol
        B = self.new_temp_symbol(blocknames[blockdata.block_mode])

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
        else:
            # The weight doesn't actually depend on anything
            # other than which quadrature rule it belongs to
            weight = self.backend.access.quadrature_weight(None, None, None, num_points)

        if blockdata.block_mode in ("full", "partial"):
            # Collect modified arguments
            ma_data = blockdata.ma_data
            rank = len(ma_data)
            assert rank == block_rank
            mts = []
            for i in range(rank):
                if blockdata.block_mode == "partial" and i == blockdata.piecewise_ma_index:
                    # Skip the piecewise index, to be applied after loop
                    continue

                mad = ma_data[i]
                td = mad.tabledata

                if td.ttype in ("quadrature", "ones"):
                    # Skip factors that contribute as 1.0
                    # FIXME: This assumes that quadrature ttype is correctly exploited in indexing!
                    continue

                if td.is_piecewise:
                    # ma_index refers to piecewise scope  # TODO: Make this an attribute of blockdata
                    scope = self.ir["piecewise_ir"]["modified_arguments"]
                else:
                    # ma_index refers to current varying scope
                    scope = self.ir["varying_irs"][num_points]["modified_arguments"]

                mts.append((i, scope[mad.ma_index], td))

            # Translate modified terminals to code
            arg_factors = []
            for i, mt, td in mts:
                access = self.backend.access(mt.terminal,
                    mt, td, num_points)
                arg_factors.append(access)

        if blockdata.block_mode == "full":
            # Multiply collected factors
            B_rhs = L.float_product([f, weight] + arg_factors)

            # Add result to block inside quadloop
            body = L.AssignAdd(B[B_indices], B_rhs)  # NB! += not =
            for i in reversed(range(block_rank)):
                if ttypes[i] != "quadrature":
                    body = L.ForRange(B_indices[i], 0, dofblock[i][1] - dofblock[i][0], body=body)
            quadparts += [body]

        elif blockdata.block_mode == "partial":
            # Multiply collected factors
            P_rhs = L.float_product([f, weight] + arg_factors)

            # Get indices and dimensions right here...  # FIXME: Validate that this transposes correctly!
            assert rank == 2
            i = blockdata.piecewise_ma_index
            not_piecewise_index = 1 - i
            P_range = dofblock[not_piecewise_index]
            P_dim = P_range[1] - P_range[0]
            P_index = self.backend.symbols.argument_loop_index(not_piecewise_index)

            # Declare P table in preparts
            P = self.new_temp_symbol("PM")
            preparts.append(L.ArrayDecl("double", P, P_dim, 0,
                                        alignas=self.ir["alignas"]))

            # Accumulate P += weight * f * args in quadrature loop
            body = L.AssignAdd(P[P_index], P_rhs)
            if ttypes[i] != "quadrature":  # FIXME: What does this mean?
                body = L.ForRange(P_index, 0, P_dim, body=body)
            quadparts.append(body)

            # Get access to element table for piecewise argument
            mt = self.ir["piecewise_ir"]["modified_arguments"][ma_data[i].ma_index]
            piecewise_argument = self.backend.access(mt.terminal,
                mt, ma_data[i].tabledata, num_points)

            # Define B = B_rhs = piecewise_argument[:] * P[:], where P[:] = sum_q weight * f * other_argument[:]
            B_rhs = piecewise_argument * P[P_index]

        if blockdata.block_mode in ("premultiplied", "preintegrated"):
            # Currently not handling facet-facet combinations for precomputed blocks
            assert self.ir["integral_type"] != "interior_facet"
            # Get the current cell entity
            entity = self.backend.symbols.entity(self.ir["entitytype"], None)

        if blockdata.block_mode == "premultiplied":
            # Declare FI in preparts
            FI = self.new_temp_symbol("FI")
            preparts.append(L.VariableDecl("double", FI, 0))

            # Accumulate FI += weight * f in quadparts
            quadparts.append(L.AssignAdd(FI, L.float_product([weight, f])))

            # Define B = B_rhs = FI * P where FI = sum_q weight*f, and P = u * v
            B_rhs = L.float_product([FI, L.Symbol(blockdata.name)[(entity,) + B_indices]])

        elif blockdata.block_mode == "preintegrated":
            # Preintegrated should never get into quadloops
            assert num_points == None

            # Define B = B_rhs = f * P where P = sum_q weight * u * v
            B_rhs = L.float_product([f, L.Symbol(blockdata.name)[(entity,) + B_indices]])

        # For full, partial, premultiplied, quadloop is non-empty
        if blockdata.block_mode in ("full", "partial", "premultiplied") and not quadparts:
            error("Not expecting empty quadparts for block mode %s." % (blockdata.block_mode,))

        # For partial, premultiplied, and preintegrated, block write occurs after quadloop
        if blockdata.block_mode != "full":
            # Write result to block
            body = L.Assign(B[B_indices], B_rhs)  # NB! = not +=
            for i in reversed(range(block_rank)):
                if ttypes[i] != "quadrature":
                    body = L.ForRange(B_indices[i], 0, dofblock[i][1] - dofblock[i][0], body=body)
            postparts.append(body)

        return B, preparts, quadparts, postparts

                # FIXME:
                # Implement this initially:
                # 1) Compute B in quadloop
                # B = 0  # Storage: num_dofs * num_dofs  // Reuse space for all factors?
                # for (q)
                #     (1a)
                #     f = factor * weight
                #     (1b)
                #     for (i)
                #         C[i] = u[i]*f
                #     (1c)
                #     for (i)
                #         for (j)
                #             B[i,j] += C[i]*v[j]
                # 2) A += B;               add block to A in finalization

                # Alternative possible optimization, needs more temporary storage:
                # 1) Compute f in quadloop    # Storage per factor: num_points
                # for (q)
                #     f1[q] = factor1 * weight[q]
                #     f#[q] = ...
                #     fn[q] = factorn * weight[q]
                # 2) Compute C in quadloop    # Storage: num_points*num_dofs
                # for (q)
                #     for (i)
                #         C[i,q] = u[i,q]*f[q]
                # 3) Compute B in quadloop    # Storage: num_dofs*num_dofs
                # B = 0
                # for (q)
                #     for (i)
                #         for (j)
                #             B[i,j] += C[i,q]*v[j,q]
                # 4) A += B;               add block to A in finalization
                # Note that storage for C and B can be reused for all terms
                # if 2-3-4 are completed for each term before starting the next.
                # Of course reusing C and B across terms where possible.


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
            A_indices = [self.backend.symbols.argument_loop_index(i)
                         for i in range(rank)]

            for dofblock, contributions in sorted(self.finalization_blocks.items()):
                # Offset indices from A indices to access dense block B
                B_indices = tuple(A_indices[i] - dofblock[i][0]
                                  for i in range(rank))

                # Sum up all blocks contributing to the same dofblock
                term = L.Sum([B[B_indices] for B in contributions])

                # Add components of all B's to A component in loop nest
                body = L.AssignAdd(A[A_indices], term)
                for i in reversed(range(rank)):
                    # TODO: need ttypes associated with this B to deal
                    # with loop dropping for quadrature elements:
                    # if ttypes[i] != "quadrature":
                    body = L.ForRange(A_indices[i], dofblock[i][0], dofblock[i][1], body=body)

                # Add this block to parts
                parts.append(body)

        return parts
