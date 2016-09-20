# -*- coding: utf-8 -*-
# Copyright (C) 2011-2015 Martin Sandve Alnæs
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

from uflacs.analysis.modified_terminals import analyse_modified_terminal, is_modified_terminal


class IntegralGenerator(object):

    def __init__(self, ir, backend):
        # Store ir
        self.ir = ir

        # Consistency check on quadrature rules
        nps1 = sorted(ir["uflacs"]["expr_irs"].keys())
        nps2 = sorted(ir["quadrature_rules"].keys())
        if nps1 != nps2:
            warning("Got different num_points for expression irs and quadrature rules:\n{0}\n{1}".format(
                nps1, nps2))

        # Compute shape of element tensor
        if self.ir["integral_type"] == "interior_facet":
            self._A_shape = [2 * n for n in self.ir["prim_idims"]]
        else:
            self._A_shape = self.ir["prim_idims"]

        #self._using_names = set()
        #self._includes = set()
        self._ufl_names = set()

        # TODO: Get self.alignas, self.padlen from ir
        sizeof_double = 8
        self.alignas = 32
        self.padlen = self.alignas // sizeof_double

        # Backend specific plugin with attributes
        # - language: for translating ufl operators to target language
        # - definitions: for defining backend specific variables
        # - access: for accessing backend specific variables
        self.backend = backend


    def generate_using_statements(self):
        #L = self.backend.language
        return []  # [L.Using(name) for name in sorted(self._using_names)]


    def get_includes(self):
        includes = set()  # self._includes)

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

        includes.update(self.backend.definitions.get_includes())

        return sorted(includes)


    def generate(self):
        """Generate entire tabulate_tensor body.

        Assumes that the code returned from here will be wrapped in a context
        that matches a suitable version of the UFC tabulate_tensor signatures.
        """
        L = self.backend.language

        parts = []
        parts += self.generate_using_statements()
        parts += self.backend.definitions.initial()
        parts += self.generate_quadrature_tables()
        parts += self.generate_element_tables()
        parts += self.generate_tensor_reset()

        # If we have integrals with different number of quadrature points,
        # we wrap each integral in a separate scope, avoiding having to
        # think about name clashes for now. This is a bit wasteful in that
        # piecewise quantities are not shared, but at least it should work.
        expr_irs = self.ir["uflacs"]["expr_irs"]
        all_num_points = sorted(expr_irs)

        # Reset variables, separate sets for quadrature loop
        self.vaccesses = { num_points: {} for num_points in all_num_points }

        for num_points in all_num_points:
            pp = self.generate_piecewise_partition(num_points)
            ql = self.generate_quadrature_loops(num_points)
            if len(all_num_points) > 1:
                # Wrapping in Scope to avoid thinking about scoping issues
                parts += [L.Scope([pp, ql])]
            else:
                parts += [pp, ql]

        parts += self.generate_finishing_statements()

        return L.StatementList(parts)


    def generate_quadrature_tables(self):
        "Generate static tables of quadrature points and weights."
        L = self.backend.language

        parts = []

        # No quadrature tables for custom (given argument) or point (evaluation in single vertex)
        if self.ir["integral_type"] in ("custom", "vertex"):
            return parts

        qrs = self.ir["quadrature_rules"]

        for num_points in sorted(qrs):
            weights = qrs[num_points][0]
            points = qrs[num_points][1]

            # Size of quadrature points depends on context, assume this is correct:
            pdim = len(points[0])

            wname = self.backend.access.weights_array_name(num_points)
            pname = self.backend.access.points_array_name(num_points)

            # Quadrature weights array
            parts += [L.ArrayDecl("static const double", wname, num_points, weights,
                                  alignas=self.alignas)]
            # Quadrature points array
            expr_ir = self.ir["uflacs"]["expr_irs"][num_points]
            if expr_ir["need_points"] and pdim > 0:
                # Flatten array: (TODO: avoid flattening here, it makes padding harder)
                points = points.reshape(product(points.shape))
                parts += [L.ArrayDecl("static const double", pname, num_points * pdim, points,
                                      alignas=self.alignas)]

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
        expr_irs = self.ir["uflacs"]["expr_irs"]
        for num_points in sorted(expr_irs):
            # Get all unique tables for this quadrature rule
            tables = expr_irs[num_points]["unique_tables"]

            # Produce code
            if len(tables):
                tmp = "Definitions of {0} tables for {1} quadrature points"
                comment = tmp.format(len(tables), num_points)
                parts += [L.Comment(comment)]
                for name in sorted(tables):
                    table = tables[name]
                    parts += [L.ArrayDecl("static const double", name, table.shape, table,
                                          alignas=self.alignas)]  # TODO: Not padding, consider when and if to do so
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
            code = tmp.format(ptrname=ptrname, size=size)
            return L.VerbatimStatement(code)

        # Compute tensor size
        A_size = product(self._A_shape)
        A = self.backend.access.element_tensor_name()

        # Stitch it together
        parts = []
        parts += [L.Comment("Reset element tensor")]
        parts += [memzero(A, A_size)]
        return parts


    def generate_quadrature_loops(self, num_points):
        "Generate all quadrature loops."
        L = self.backend.language
        parts = []

        body = self.generate_quadrature_body(num_points)
        iq = self.backend.access.quadrature_loop_index()

        if num_points == 1:
            # Wrapping body in Scope to avoid thinking about scoping issues
            # TODO: Specialize generated code with iq=0 instead of defining iq here.
            parts += [L.Comment("Only 1 quadrature point, no loop"),
                      L.VariableDecl("const int", iq, 0),
                      L.Scope(body)]

        else:
            parts += [L.ForRange(iq, 0, num_points, body=body)]
        return parts


    def generate_quadrature_body(self, num_points):
        """
        """
        parts = []
        L = self.backend.language
        parts += self.generate_varying_partition(num_points)
        if parts:
            parts = [L.Comment("Quadrature loop body setup (num_points={0})".format(num_points))] + parts

        # Compute single argument partitions outside of the dofblock loops
        for iarg in range(self.ir["rank"]):
            for dofrange in []:  # TODO: Move f*arg0 out here
                parts += self.generate_argument_partition(num_points, iarg, dofrange)

        # Nested argument loops and accumulation into element tensor
        parts += self.generate_quadrature_body_dofblocks(num_points)

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

        expr_ir = self.ir["uflacs"]["expr_irs"][num_points]
        # tuple(modified_argument_indices) -> code_index
        AF = expr_ir["argument_factorization"]

        # modified_argument_index -> (tablename, dofbegin, dofend)
        MATR = expr_ir["modified_argument_table_ranges"]

        # Find dofranges at this loop level iarg starting with outer_dofblock
        dofranges = set()
        for mas in AF:
            mas_full_dofblock = tuple(MATR[j][1:3] for j in mas)
            if tuple(mas_full_dofblock[:iarg]) == tuple(outer_dofblock):
                dofrange = mas_full_dofblock[iarg]
                assert dofrange[0] != dofrange[1]
                dofranges.add(dofrange)
        dofranges = sorted(dofranges)

        # Build loops for each dofrange
        for dofrange in dofranges:
            dofblock = outer_dofblock + (dofrange,)

            # Generate nested inner loops (only triggers for forms with two or more arguments
            body = self.generate_quadrature_body_dofblocks(num_points, dofblock)

            # Wrap setup, subloops, and accumulation in a loop for this level
            idof = self.backend.access.argument_loop_index(iarg)
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
                vaccess = self.backend.access(mt.terminal, mt, table_ranges[i], num_points)
                vdef = self.backend.definitions(mt.terminal, mt, table_ranges[i], vaccess)

                # Store definitions of terminals in list
                if vdef is not None:
                    definitions.append(vdef)
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

                # Currently instead creating a new intermediate for each subexpression except boolean conditions
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

        # Join terminal computation, array of intermediate expressions, and intermediate computations
        parts = [definitions]
        if intermediates:
            parts += [L.ArrayDecl("double", symbol, len(intermediates), alignas=self.alignas)]
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

                if isinstance(mt.terminal, ConstantValue):
                    # Format literal value for the chosen language
                    vaccess = modified_literal_to_ast_node[tc](mt)  # XXX FIXME: Implement this mapping
                    vdef = None
                else:
                    # Backend specific modified terminal formatting
                    vaccess = self.backend.access(mt.terminal, mt, table_ranges[i], num_points)
                    vdef = self.backend.definitions(mt.terminal, mt, table_ranges[i], vaccess)

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

        # Join terminal computation, array of intermediate expressions, and intermediate computations
        parts = [definitions]
        if intermediates:
            parts += [L.ArrayDecl("double", symbol, len(intermediates), alignas=self.alignas)]
            parts += intermediates
        return parts


    def generate_piecewise_partition(self, num_points):
        """Generate statements prior to the quadrature loop.

        This mostly includes computations involving piecewise constant geometry and coefficients.
        """
        L = self.backend.language
        expr_ir = self.ir["uflacs"]["expr_irs"][num_points]
        arraysymbol = L.Symbol("sp{0}".format(num_points))
        parts = self.generate_partition(arraysymbol,
                                        expr_ir["V"],
                                        expr_ir["piecewise"],
                                        expr_ir["table_ranges"],
                                        num_points)
        if parts:
            parts.insert(0, L.Comment("Section for piecewise constant computations"))
        return parts


    def generate_varying_partition(self, num_points):
        L = self.backend.language
        expr_ir = self.ir["uflacs"]["expr_irs"][num_points]
        arraysymbol = L.Symbol("sv{0}".format(num_points))
        parts = self.generate_partition(arraysymbol,
                                        expr_ir["V"],
                                        expr_ir["varying"],
                                        expr_ir["table_ranges"],
                                        num_points)
        if parts:
            parts.insert(0, L.Comment("Section for geometrically varying computations"))
        return parts


    def generate_argument_partition(self, num_points, iarg, dofrange):
        """Generate code for the partition corresponding to arguments 0..iarg within given dofblock."""
        parts = []
        # TODO: What do we want to do here? Define!
        # Should this be a single loop over i0, i1 separately
        # outside of the double loop over (i0,i1)?
        return parts


    def generate_integrand_accumulation(self, num_points, dofblock):
        parts = []
        L = self.backend.language

        expr_ir = self.ir["uflacs"]["expr_irs"][num_points]
        AF = expr_ir["argument_factorization"]
        V = expr_ir["V"]
        MATR = expr_ir["modified_argument_table_ranges"]
        MA = expr_ir["modified_arguments"]

        idofs = [self.backend.access.argument_loop_index(i) for i in range(self.ir["rank"])]

        # Find the blocks to build: (TODO: This is rather awkward,
        # having to rediscover these relations here)
        arguments_and_factors = sorted(expr_ir["argument_factorization"].items(),
                                       key=lambda x: x[0])
        for args, factor_index in arguments_and_factors:
            if not all(tuple(dofblock[iarg]) == tuple(MATR[ma][1:3])
                       for iarg, ma in enumerate(args)):
                continue

            factors = []

            # Get factor expression
            v = V[factor_index]
            if v._ufl_is_literal_ and float(v) == 1.0:
                # TODO: Nicer way to check for f=1?
                pass
            else:
                fexpr = self.vaccesses[num_points][v]
                factors.append(fexpr)

            # Get table names
            argfactors = []
            for i, ma in enumerate(args):
                access = self.backend.access(MA[ma].terminal, MA[ma], MATR[ma], num_points)
                argfactors += [access]

            factors.extend(argfactors)

            # Format index access to A
            A_access = self.backend.access.element_tensor_entry(idofs, self._A_shape)

            # Emit assignment
            parts += [L.AssignAdd(A_access, L.Product(factors))]

        return parts


    def generate_finishing_statements(self):
        """Generate finishing statements.

        This includes assigning to output array if there is no integration.
        """
        parts = []

        if not self.ir["quadrature_rules"]:  # Rather check ir["integral_type"]?
            # TODO: Implement for expression support
            error("Expression generation not implemented yet.")
            # TODO: If no integration, assuming we generate an expression, and assign results here
            # Corresponding code from compiler.py:
            # assign_to_variables = tfmt.output_variable_names(len(final_variable_names))
            # parts += list(format_assignments(zip(assign_to_variables, final_variable_names)))

        return parts
