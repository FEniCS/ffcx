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

from ufl import product
from ufl.classes import ConstantValue, Condition

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


    def generate(self):
        """Generate entire tabulate_tensor body.

        Assumes that the code returned from here will be wrapped in a context
        that matches a suitable version of the UFC tabulate_tensor signatures.
        """
        L = self.backend.language

        parts = []
        parts += self.generate_quadrature_tables()
        parts += self.generate_element_tables()
        parts += self.generate_tensor_reset()

        # If we have integrals with different number of quadrature points,
        # we wrap each integral in a separate scope, avoiding having to
        # think about name clashes for now. This is a bit wasteful in that
        # piecewise quantities are not shared, but at least it should work.
        expr_irs = self.ir["expr_irs"]
        all_num_points = sorted(expr_irs)

        # Reset variables, separate sets for quadrature loop
        self.vaccesses = { num_points: {} for num_points in all_num_points }

        for num_points in all_num_points:
            pp = self.generate_piecewise_partition(num_points)
            ql = self.generate_quadrature_loops(num_points)
            body = pp + ql
            if len(all_num_points) > 1:
                # Wrapping in Scope to avoid thinking about scoping issues
                parts.append(L.Scope(body))
            else:
                parts.extend(body)

        parts += self.generate_finishing_statements()

        return L.StatementList(parts)


    def generate_quadrature_tables(self):
        "Generate static tables of quadrature points and weights."
        L = self.backend.language

        parts = []

        # No quadrature tables for custom (given argument)
        # or point (evaluation in single vertex)
        skip = ("custom", "cutcell", "interface", "overlap", "vertex")
        if self.ir["integral_type"] in skip:
            return parts

        # Loop over quadrature rules
        qrs = self.ir["quadrature_rules"]
        for num_points in sorted(qrs):
            points, weights = qrs[num_points]
            assert num_points == len(weights)
            expr_ir = self.ir["expr_irs"][num_points]

            # Generate quadrature weights array
            if expr_ir["need_weights"]:
                wsym = self.backend.symbols.weights_array(num_points)
                parts += [L.ArrayDecl("static const double", wsym, num_points, weights,
                                      alignas=self.ir["alignas"])]

            # Size of quadrature points depends on context, assume this is correct:
            pdim = len(points[0])
            assert points.shape[0] == num_points
            assert pdim == points.shape[1]
            #import IPython; IPython.embed()

            # Generate quadrature points array
            if pdim and expr_ir["need_points"]:
                # Flatten array: (TODO: avoid flattening here, it makes padding harder)
                flattened_points = points.reshape(product(points.shape))
                psym = self.backend.symbols.points_array(num_points)
                parts += [L.ArrayDecl("static const double", psym, num_points * pdim,
                                      flattened_points, alignas=self.ir["alignas"])]

        # Add leading comment if there are any tables
        if parts:
            header = [L.Comment("Section for quadrature weights and points")]
            parts = header + parts
        return parts


    def generate_element_tables(self):
        """Generate static tables with precomputed element basis
        function values in quadrature points."""
        L = self.backend.language
        parts = []
        expr_irs = self.ir["expr_irs"]

        for num_points in sorted(expr_irs):
            # Get all unique tables for this quadrature rule
            tables = expr_irs[num_points]["unique_tables"]
            if tables:
                tmp = "Definitions of {0} tables for {1} quadrature points"
                parts += [L.Comment(tmp.format(len(tables), num_points))]
                for name in sorted(tables):
                    # TODO: table here can actually have only 1 point,
                    # regroup or at least fix generated comment
                    table = tables[name]
                    # TODO: Not padding, consider when and if to do so
                    parts += [L.ArrayDecl("static const double", name, table.shape, table,
                                          alignas=self.ir["alignas"])]
        # Add leading comment if there are any tables
        if parts:
            header = [L.Comment("Section for precomputed element basis function values"),
                      L.Comment("Table dimensions: num_entities, num_points, num_dofs")]
            parts = header + parts
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


    def generate_quadrature_loops(self, num_points):
        "Generate all quadrature loops."
        L = self.backend.language
        body = []

        # Generate unstructured varying partition
        body += self.generate_varying_partition(num_points)
        if body:
            comment = [L.Comment("Quadrature loop body setup (num_points={0})".format(num_points))]
            body = comment + body

        """Previous data flow:
        generate_piecewise_partition
        generate_quadrature_loops
         generate_quadrature_body
          generate_varying_partition
           generate_argument_partition
            generate_quadrature_body_dofblocks
             generate_quadrature_body_dofblocks
              generate_integrand_accumulation
        """

        if 0:
            # Nested argument loops and accumulation into element tensor
            body += self.generate_quadrature_body_dofblocks(num_points)
        else:
            # FIXME: Use this instead:
            body += self.generate_quadrature_body_dofblocks2(num_points)

        # Wrap body in loop or scope
        if not body:
            # Could happen for integral with everything zero and optimized away
            parts = []
        elif num_points == 1:
            # For now wrapping body in Scope to avoid thinking about scoping issues
            parts = [L.Comment("Only 1 quadrature point, no loop"),
                     L.Scope(body)]
        else:
            # Regular case: define quadrature loop
            iq = self.backend.symbols.quadrature_loop_index(num_points)
            np = self.backend.symbols.num_quadrature_points(num_points)
            parts = [L.ForRange(iq, 0, np, body=body)]

        return parts


    def generate_quadrature_body_dofblocks2(self, num_points):
        L = self.backend.language

        # TODO: Add partial blocks (T[i0] = factor_index * arg0;)

        # TODO: Move piecewise blocks outside quadrature loop
        # (Can only do this by removing weight from factor,
        # and using that piecewise f*u*v gives that
        # sum_q weight[q]*f*u*v == f*u*v*(sum_q weight[q]) )

        # Get representation details
        expr_ir = self.ir["expr_irs"][num_points]
        V = expr_ir["V"]
        modified_arguments = expr_ir["modified_arguments"]
        block_contributions = expr_ir["block_contributions"]

        vaccesses = self.vaccesses[num_points]
        A = self.backend.symbols.element_tensor()

        parts = []
        for partition in ["piecewise", "varying"]:
            partition_parts = []
            for dofblock, contributions in sorted(block_contributions[partition].items()):
                for data in contributions:
                    (ma_indices, factor_index, table_ranges, unames, ttypes) = data

                    # Add code in layers starting with innermost A[...] += product(factors)
                    rank = len(unames)
                    factors = []

                    # Get factor expression
                    v = V[factor_index]
                    if not (v._ufl_is_literal_ and float(v) == 1.0):
                        factors.append(vaccesses[v])

                    # Get loop counter symbols to access A with
                    A_indices = []
                    for i in range(rank):
                        if ttypes[i] == "quadrature":
                            # Used to index A like A[iq*num_dofs + iq]
                            ia = self.backend.symbols.quadrature_loop_index(num_points)
                        else:
                            # Regular dof index
                            ia = self.backend.symbols.argument_loop_index(i)
                        A_indices.append(ia)

                    # Add table access to factors, unless it's always 1.0
                    for i in range(rank):
                        tt = ttypes[i]
                        assert tt not in ("zeros",)
                        if tt not in ("quadrature", "ones"):
                            ma = ma_indices[i]
                            access = self.backend.access(
                                modified_arguments[ma].terminal,
                                modified_arguments[ma],
                                table_ranges[i],
                                num_points)
                            factors.append(access)

                    # Special case where all factors are 1.0 and dropped
                    if factors:
                        term = L.Product(factors)
                    else:
                        term = L.LiteralFloat(1.0)

                    # Format flattened index expression to access A
                    flat_index = L.flattened_indices(A_indices, self.ir["tensor_shape"])
                    body = [L.AssignAdd(A[flat_index], term)]

                    # Wrap accumulation in loop nest
                    #for i in range(rank):
                    for i in range(rank-1, -1, -1):
                        if ttypes[i] != "quadrature":
                            dofrange = dofblock[i]
                            body = [L.ForRange(A_indices[i], dofrange[0], dofrange[1], body=body)]

                    # Add this block to parts
                    partition_parts.append(body)
            if partition_parts:
                parts += [L.Comment("Blocks of %s contributions" % (partition,))]
                parts.extend(partition_parts)
        return parts


    def generate_quadrature_body_dofblocks(self, num_points, outer_dofblock=()):
        parts = []
        L = self.backend.language

        # The loop level iarg here equals the argument count (in renumbered >= 0 format)
        iarg = len(outer_dofblock)
        if iarg == self.ir["rank"]:
            # At the innermost argument loop level we accumulate into the element tensor
            parts += [self.generate_integrand_accumulation(num_points, outer_dofblock)]
            return parts
        assert iarg < self.ir["rank"]

        expr_ir = self.ir["expr_irs"][num_points]

        # tuple(modified_argument_indices) -> code_index
        argument_factorization = expr_ir["argument_factorization"]

        # modified_argument_index -> (tablename, dofbegin, dofend)
        modified_argument_table_ranges = expr_ir["modified_argument_table_ranges"]

        # Find dofranges at this loop level iarg starting with outer_dofblock
        dofranges = set()
        for ma_indices in argument_factorization:
            ma_table_ranges = tuple(modified_argument_table_ranges[j]
                                    for j in ma_indices)
            full_dofblock = tuple(tr[1:3] for tr in ma_table_ranges)

            # FIXME: Figure out exactly why iarg < len(full_dofblock) is needed
            if iarg < len(full_dofblock) and tuple(full_dofblock[:iarg]) == tuple(outer_dofblock):
                dofrange = full_dofblock[iarg]
                assert dofrange[0] != dofrange[1]
                dofranges.add(dofrange)
        dofranges = sorted(dofranges)

        # Build loops for each dofrange
        for dofrange in dofranges:
            dofblock = outer_dofblock + (dofrange,)

            # Generate nested inner loops (only triggers for forms with two or more arguments
            body = self.generate_quadrature_body_dofblocks(num_points, dofblock)

            # Wrap setup, subloops, and accumulation in a loop for this level
            ttype = "not-quadrature"  # FIXME
            if ttype == "quadrature":
                # FIXME: Possibly need to extract the quadrature dofranges at a higher level
                # For quadrature element, drop the loop
                parts += [body]
            else:
                idof = self.backend.symbols.argument_loop_index(iarg)
                parts += [L.ForRange(idof, dofrange[0], dofrange[1], body=body)]
        return parts


    def generate_partition(self, symbol, V, partition, table_ranges, num_points):
        L = self.backend.language

        definitions = []
        intermediates = []

        vaccesses = self.vaccesses[num_points]

        partition_indices = [i for i, p in enumerate(partition) if p]

        for i in partition_indices:
            v = V[i]

            if is_modified_terminal(v):
                mt = analyse_modified_terminal(v)

                # Backend specific modified terminal translation
                vaccess = self.backend.access(mt.terminal,
                    mt, table_ranges[i], num_points)
                vdef = self.backend.definitions(mt.terminal,
                    mt, table_ranges[i], num_points, vaccess)

                # Store definitions of terminals in list
                assert isinstance(vdef, list)
                definitions.extend(vdef)
            else:
                # Get previously visited operands (TODO: use edges of V instead of ufl_operands?)
                vops = [vaccesses[op] for op in v.ufl_operands]

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
            vaccesses[v] = vaccess

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


    # TODO: Rather take list of vertices, not markers
    # XXX FIXME: Fix up this function and use it instead!
    def alternative_generate_partition(self, symbol, C, MT, partition, table_ranges, num_points):
        L = self.backend.language

        definitions = []
        intermediates = []

        # XXX FIXME: create these!
        # C = input CRSArray representation of expression DAG
        # MT = input list/dict of modified terminals

        self.ast_variables = [None]*len(C) # FIXME: Create outside

        # TODO: Get this as input instead of partition?
        partition_indices = [i for i, p in enumerate(partition) if p]
        for i in partition_indices:
            row = C[i] # XXX FIXME: Get this as input
            if len(row) == 1:
                # Modified terminal
                t, = row
                mt = MT[t] # XXX FIXME: Get this as input
                tc = mt[0]


                if isinstance(mt.terminal, ConstantValue):
                    # Format literal value for the chosen language
                    modified_literal_to_ast_node = []  # silence flake8
                    # XXX FIXME: Implement this mapping:
                    vaccess = modified_literal_to_ast_node[tc](mt)
                    vdef = None
                else:
                    # Backend specific modified terminal formatting
                    vaccess = self.backend.access(mt.terminal,
                        mt, table_ranges[i], num_points)
                    vdef = self.backend.definitions(mt.terminal,
                        mt, table_ranges[i], num_points, vaccess)

                # Store definitions of terminals in list
                if vdef is not None:
                    definitions.append(vdef)

            else:
                # Application of operator with typecode tc to operands with indices ops
                tc = mt[0]
                ops = mt[1:]

                # Get operand AST nodes
                opsaccess = [self.ast_variables[k] for k in ops]

                # Generate expression for this operator application
                typecode2astnode = []  # silence flake8
                vexpr = typecode2astnode[tc](opsaccess) # XXX FIXME: Implement this mapping

                store_this_in_variable = True # TODO: Don't store all subexpressions
                if store_this_in_variable:
                    # Record assignment of vexpr to intermediate variable
                    j = len(intermediates)
                    vaccess = symbol[j]
                    intermediates.append(L.Assign(vaccess, vexpr))
                else:
                    # Access the inlined expression
                    vaccess = vexpr

            # Store access string, either a variable symbol or an inlined expression
            self.ast_variables[i] = vaccess

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


    def generate_piecewise_partition(self, num_points):
        """Generate statements prior to the quadrature loop.

        This mostly includes computations involving piecewise constant geometry and coefficients.
        """
        L = self.backend.language
        expr_ir = self.ir["expr_irs"][num_points]
        arraysymbol = L.Symbol("sp{0}".format(num_points))
        parts = self.generate_partition(arraysymbol,
                                        expr_ir["V"],
                                        expr_ir["piecewise"],
                                        expr_ir["table_ranges"],
                                        num_points)
        if parts:
            comment = [L.Comment("Section for piecewise constant computations")]
            parts = comment + parts
        return parts


    def generate_varying_partition(self, num_points):
        L = self.backend.language
        expr_ir = self.ir["expr_irs"][num_points]
        arraysymbol = L.Symbol("sv{0}".format(num_points))
        parts = self.generate_partition(arraysymbol,
                                        expr_ir["V"],
                                        expr_ir["varying"],
                                        expr_ir["table_ranges"],
                                        num_points)
        if parts:
            comment = [L.Comment("Section for geometrically varying computations")]
            parts = comment + parts
        return parts


    def generate_integrand_accumulation(self, num_points, dofblock):
        parts = []
        L = self.backend.language

        vaccesses = self.vaccesses[num_points]

        # Get representation details
        expr_ir = self.ir["expr_irs"][num_points]
        argument_factorization = expr_ir["argument_factorization"]
        V = expr_ir["V"]
        modified_argument_table_ranges = expr_ir["modified_argument_table_ranges"]
        modified_arguments = expr_ir["modified_arguments"]

        # Get some symbols
        A = self.backend.symbols.element_tensor()
        idofs = []
        for i in range(self.ir["rank"]):
            # FIXME: If ttype == "quadrature", drop loop and set argument access == 1.0
            ttype = "not-quadrature"  # FIXME
            if ttype == "quadrature":
                # Used to index A like A[iq*num_dofs + iq]
                ia = self.backend.symbols.quadrature_loop_index(num_points)
            else:
                ia = self.backend.symbols.argument_loop_index(i)
            idofs.append(ia)

        # Find the blocks to build: (TODO: This is rather awkward,
        # having to rediscover these relations here)
        arguments_and_factors = sorted(argument_factorization.items(),
                                       key=lambda x: x[0])
        for args, factor_index in arguments_and_factors:
            if not all(tuple(dofblock[iarg])
                       == tuple(modified_argument_table_ranges[ma][1:3])
                       for iarg, ma in enumerate(args)):
                continue

            factors = []

            # Get factor expression
            v = V[factor_index]
            # TODO: Nicer way to check for f=1?
            if not (v._ufl_is_literal_ and float(v) == 1.0):
                factors.append(vaccesses[v])

            # Get table names
            for i, ma in enumerate(args):
                # FIXME: If ttype == "quadrature", drop loop and set argument access == 1.0
                ttype = "not-quadrature"  # FIXME
                if ttype not in ("quadrature", "ones"):
                    access = self.backend.access(
                        modified_arguments[ma].terminal,
                        modified_arguments[ma],
                        modified_argument_table_ranges[ma],
                        num_points)
                    factors.append(access)

            # Special case where all factors are 1.0 and dropped
            if factors:
                term = L.Product(factors)
            else:
                term = L.LiteralFloat(1.0)

            # Format flattened index expression to access A
            flat_index = L.flattened_indices(idofs, self.ir["tensor_shape"])

            # Emit assignment
            parts += [L.AssignAdd(A[flat_index], term)]

        return parts


    def generate_finishing_statements(self):
        """Generate finishing statements.

        This includes assigning to output array if there is no integration.
        """
        parts = []

        if self.ir["integral_type"] == "expression":
            error("Expression generation not implemented yet.")
            # TODO: If no integration, assuming we generate an expression, and assign results here
            # Corresponding code from compiler.py:
            # assign_to_variables = tfmt.output_variable_names(len(final_variable_names))
            # parts += list(format_assignments(zip(assign_to_variables, final_variable_names)))

        return parts
