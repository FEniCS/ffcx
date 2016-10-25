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

        # Reset variables, separate sets for quadrature loop
        self.vaccesses = { num_points: {} for num_points in self.ir["all_num_points"] }
        self.vaccesses[None] = {}


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


    def get_vaccess(self, v, num_points):
        # Lookup v in 'num_points' quadloop scope
        f = self.vaccesses[num_points].get(v)
        if f is None:
            # Missed quadloop scope lookup, lookup in piecewise scope
            # (this should exist now)
            f = self.vaccesses[None][v]
        else:
            #assert num_points is None or v not in self.vaccesses[None]
            v in self.vaccesses[None]
        return f


    def generate(self):
        """Generate entire tabulate_tensor body.

        Assumes that the code returned from here will be wrapped in a context
        that matches a suitable version of the UFC tabulate_tensor signatures.
        """
        L = self.backend.language

        assert not any(d for d in self.vaccesses.values())

        parts = []
        parts += self.generate_quadrature_tables()
        parts += self.generate_element_tables()
        parts += self.generate_tensor_reset()
        parts += self.generate_piecewise_partition()
        parts += self.generate_preintegrated_dofblocks()
        parts += self.generate_piecewise_dofblock_partition()

        # If we have integrals with different number of quadrature points,
        # we wrap each integral in a separate scope, avoiding having to
        # think about name clashes for now. This is a bit wasteful in that
        # piecewise quantities are not shared, but at least it should work.
        all_num_points = self.ir["all_num_points"]
        for num_points in all_num_points:
            body = self.generate_quadrature_loops(num_points)

            # If there are multiple quadrature rules here, just wrapping
            # in Scope to avoid thinking about scoping issues for now.
            # A better handling of multiple rules would be nice,
            # in particular 
            if len(all_num_points) > 1:
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
            "Precomputed values of basis functions", # and partial integrals",
            "FE* dimensions: [entities][points][dofs]",
            # TODO: When adding preintegrated tables extend the comments here:
            #"B* dimensions: [entities][dofs]", 
            #"C* dimensions: [entities][dofs][dofs]",
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


    def generate_quadrature_loops(self, num_points):
        "Generate all quadrature loops."
        L = self.backend.language
        body = []

        # Generate unstructured varying partition
        body += self.generate_varying_partition(num_points)
        body = L.commented_code_list(body,
            "Quadrature loop body setup (num_points={0})".format(num_points))

        body += self.generate_varying_dofblock_partition(num_points)

        # Wrap body in loop or scope
        if not body:
            # Could happen for integral with everything zero and optimized away
            parts = []
        elif num_points == 1:
            # For now wrapping body in Scope to avoid thinking about scoping issues
            parts = L.commented_code_list(L.Scope(body), "Only 1 quadrature point, no loop")
        else:
            # Regular case: define quadrature loop
            iq = self.backend.symbols.quadrature_loop_index(num_points)
            np = self.backend.symbols.num_quadrature_points(num_points)
            parts = [L.ForRange(iq, 0, np, body=body)]

        return parts


    def generate_preintegrated_dofblocks(self):
        L = self.backend.language
        A = self.backend.symbols.element_tensor()
        parts = []

        pir = self.ir["piecewise_ir"]
        #PB = pir["preintegrated_blocks"]
        PC = pir["preintegrated_contributions"]
        for dofblock, contributions in sorted(PC.items()):
            rank = len(dofblock)
            for data in contributions:
                assert data.block_mode == "preintegrate"

                v = pir["V"][data.factor_index]
                f = self.get_vaccess(v, None)

                # Get loop counter symbols to access A with
                A_indices = []
                for i in range(rank):
                    ia = self.backend.symbols.argument_loop_index(i)
                    A_indices.append(ia)

                # Offset A indices to define P indices
                restriction = None
                assert self.ir["integral_type"] != "interior_facet"
                entity = self.backend.symbols.entity(self.ir["entitytype"], restriction)
                P_indices = (entity,) + tuple(A_indices[i] - dofblock[i][0] for i in range(rank))

                # The preintegrated term scaled with piecewise factor
                P = L.Symbol(data.pname)
                term = f * P[P_indices]

                # Format flattened index expression to access A
                flat_index = L.flattened_indices(A_indices, self.ir["tensor_shape"])
                body = L.AssignAdd(A[flat_index], term)

                # Wrap accumulation in loop nest
                for i in range(rank-1, -1, -1):
                    dofrange = dofblock[i]
                    body = L.ForRange(A_indices[i], dofrange[0], dofrange[1], body=body)

                # Add this block to parts
                parts.append(body)

        # FIXME: Generate code for preintegrated_contributions
        # 1) P = weight*u*v;  preintegrate block here
        # 2) B = f*P;         scale block after quadloop
        # 3) A[dofblock] += B[:];   add block to A in finalization

        parts = L.commented_code_list(parts,
            "Preintegrated dofblocks")
        return parts


    def generate_piecewise_dofblock_partition(self):
        return self.generate_dofblock_partition(None)


    def generate_varying_dofblock_partition(self, num_points):
        return self.generate_dofblock_partition(num_points)


    def generate_dofblock_partition(self, num_points):
        L = self.backend.language
        A = self.backend.symbols.element_tensor()

        # TODO: Add partial blocks (T[i0] = factor_index * arg0;)

        # TODO: Move piecewise blocks outside quadrature loop
        # (Can only do this by removing weight from factor,
        # and using that piecewise f*u*v gives that
        # sum_q weight[q]*f*u*v == f*u*v*(sum_q weight[q]) )

        pir = self.ir["piecewise_ir"]
        if num_points is None:  # NB! meaning piecewise partition, not custom integral
            block_contributions = pir["block_contributions"]
        else:
            vir = self.ir["varying_irs"].get(num_points)
            block_contributions = vir["block_contributions"]

        parts = []
        for dofblock, contributions in sorted(block_contributions.items()):
            for blockdata in contributions:
                # TODO: Generate different code for functional, partial, and runtime
                assert blockdata.block_mode in ("functional", "partial", "runtime")
                assert blockdata.block_mode not in ("preintegrate",)
                ma_data = blockdata.ma_data
                rank = len(ma_data)

                # Add code in layers starting with innermost A[...] += product(factors)
                factors = []

                # Get factor expression
                if blockdata.factor_is_piecewise:
                    v = pir["V"][blockdata.factor_index]
                else:
                    v = vir["V"][blockdata.factor_index]
                if not (v._ufl_is_literal_ and float(v) == 1.0):
                    factors.append(self.get_vaccess(v, num_points))

                # Get loop counter symbols to access A with
                A_indices = []
                for i in range(rank):
                    if ma_data[i].tabledata.ttype == "quadrature":
                        # Used to index A like A[iq*num_dofs + iq]
                        ia = self.backend.symbols.quadrature_loop_index(num_points)
                    else:
                        # Regular dof index
                        ia = self.backend.symbols.argument_loop_index(i)
                    A_indices.append(ia)

                # Add table access to factors, unless it's always 1.0
                for i in range(rank):
                    ttype = ma_data[i].tabledata.ttype
                    # 1.0 factors do not contribute
                    if ttype in ("quadrature", "ones"):
                        pass
                    elif ttype == "zeros":
                        # Zeros should be optimized away
                        error("Not expecting zero arguments to be left in dofblock generation.")
                    else:
                        if ttype in ("piecewise", "fixed"):
                            mt = pir["modified_arguments"][ma_data[i].ma_index]
                        elif ttype in ("uniform", "varying"):
                            mt = vir["modified_arguments"][ma_data[i].ma_index]
                        else:
                            error("Not expecting table type %s in dofblock generation." % (ttype,))
                        access = self.backend.access(mt.terminal,
                            mt, ma_data[i].tabledata, num_points)
                        factors.append(access)

                # Special case where all factors are 1.0 and dropped
                if factors:
                    term = L.Product(factors)
                else:
                    term = L.LiteralFloat(1.0)

                # Format flattened index expression to access A
                flat_index = L.flattened_indices(A_indices, self.ir["tensor_shape"])
                body = L.AssignAdd(A[flat_index], term)

                # Wrap accumulation in loop nest
                #for i in range(rank):
                for i in range(rank-1, -1, -1):
                    if ma_data[i].tabledata.ttype != "quadrature":
                        dofrange = dofblock[i]
                        body = L.ForRange(A_indices[i], dofrange[0], dofrange[1], body=body)

                # Add this block to parts
                parts.append(body)

        return parts


    def generate_partition(self, symbol, V, V_active, V_table_data, num_points):
        L = self.backend.language

        definitions = []
        intermediates = []

        active_indices = [i for i, p in enumerate(V_active) if p]

        for i in active_indices:
            v = V[i]

            if is_modified_terminal(v):
                mt = analyse_modified_terminal(v)

                tabledata = V_table_data[i]

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
                vops = [self.get_vaccess(op, num_points) for op in v.ufl_operands]
                
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
            self.vaccesses[num_points][v] = vaccess

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


    def generate_piecewise_partition(self):
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


    def generate_varying_partition(self, num_points):
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
