# -*- coding: utf-8 -*-
# Copyright (C) 2013-2017 Martin Sandve Alnæs
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Controlling algorithm for building the tabulate_tensor source structure from factorized representation."""

import itertools
import functools
import logging
from copy import copy
from collections import defaultdict

from ffc import FFCError
from ffc.uflacs.build_uflacs_ir import get_common_block_data
from ffc.uflacs.elementtables import piecewise_ttypes
from ffc.uflacs.language.cnodes import pad_dim, pad_innermost_dim
from ufl import product
from ufl.classes import Condition
from ufl.measure import custom_integral_types, point_integral_types

logger = logging.getLogger(__name__)


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
        """Return list of include statements needed to support generated code."""
        includes = set()

        # Get std::fill used by MemZero
        # includes.add("#include <algorithm>")

        # For controlling floating point environment
        # includes.add("#include <cfenv>")

        # For intel intrinsics and controlling floating point environment
        # includes.add("#include <xmmintrin.h>")

        cmath_names = set((
            "abs",
            "sign",
            "pow",
            "sqrt",
            "exp",
            "ln",
            "cos",
            "sin",
            "tan",
            "acos",
            "asin",
            "atan",
            "atan_2",
            "cosh",
            "sinh",
            "tanh",
            "acosh",
            "asinh",
            "atanh",
            "erf",
            "erfc",
        ))

        boost_math_names = set((
            "bessel_j",
            "bessel_y",
            "bessel_i",
            "bessel_k",
        ))

        # Only return the necessary headers
        if cmath_names & self._ufl_names:
            includes.add("#include <math.h>")

        includes.add("#include <stdalign.h>")

        if boost_math_names & self._ufl_names:
            includes.add("#include <boost/math/special_functions.hpp>")

        return sorted(includes)

    def init_scopes(self):
        """Initialize variable scope dicts."""
        # Reset variables, separate sets for quadrature loop
        self.scopes = {num_points: {} for num_points in self.ir["all_num_points"]}
        self.scopes[None] = {}

    def set_var(self, num_points, v, vaccess):
        """Set a new variable in variable scope dicts.

        Scope is determined by num_points which identifies the
        quadrature loop scope or None if outside quadrature loops.

        v is the ufl expression and vaccess is the CNodes
        expression to access the value in the code.

        """
        self.scopes[num_points][v] = vaccess

    def has_var(self, num_points, v):
        """Check if variable exists in variable scope dicts.

        Return True if ufl expression v exists in the num_points scope.

        NB! Does not fall back to piecewise scope.
        """
        return v in self.scopes[num_points]

    def get_var(self, num_points, v):
        """Lookup ufl expression v in variable scope dicts.

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
        """Create a new code symbol named basename + running counter."""
        L = self.backend.language
        name = "%s%d" % (basename, self.symbol_counters[basename])
        self.symbol_counters[basename] += 1
        return L.Symbol(name)

    def get_temp_symbol(self, tempname, key):
        key = (tempname, ) + key
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

        # Optionally perform cross-element vectorization of the generated code
        cross_element_width = self.ir["params"]["cross_element_width"]
        if cross_element_width > 0:
            # Optionally use gcc vector extensions
            if self.ir["params"]["enable_cross_element_gcc_ext"]:
                vectorized = self.vectorize_with_gcc_exts(parts, vec_length=cross_element_width)
            else:
                vectorized = self.vectorize_with_loops(parts, vec_length=cross_element_width)

            parts = vectorized

        return L.StatementList(parts)

    def vectorize_with_gcc_exts(self, statements, vec_length=4, alignment=32):
        """
        Converts a list of tabulate_tensor CNodes statements into a cross-element vectorized version using
        GCC's vector extensions.

        :param statements: The list of statements that should be transformed
        :param vec_length: The number of elements to perform cross element vectorization over
        :param alignment: Alignment in bytes used for arrays whose rank is increased for vectorization
        :return: A list of the transformed statements
        """

        L = self.backend.language

        ctx = {
            "vectorized_inputs": ["A", "w", "coordinate_dofs"],
            "vectorized_intermediates": [],
        }

        base_type = "double"
        vector_type = "double{}".format(str(vec_length))

        # Symbol used as index in for-loops over the elements
        i_simd = L.Symbol("i_elem")

        def simd_loop(body):
            """Returns a ForRange that is used to perform a cross element loop."""
            return L.ForRange(i_simd, 0, vec_length, body)

        def was_vectorized(expr) -> bool:
            """Returns whether an array's/scalar's rank was increased for vectorization."""

            if isinstance(expr, L.ArrayAccess):
                symbol_name = expr.array.name
            elif isinstance(expr, L.Symbol):
                symbol_name = expr.name
            else:
                raise RuntimeError("Unsupported expression")

            return any(
                (symbol_name in arr) for arr in [ctx["vectorized_intermediates"], ctx["vectorized_inputs"]])

        def child_nodes(expr: L.CExpr):
            """Returns a list of all sub-expression nodes of a CNodes expression."""

            if isinstance(expr, L.CExprTerminal):
                return []
            elif isinstance(expr, L.UnaryOp):
                return [expr.arg]
            elif isinstance(expr, L.BinOp):
                return [expr.lhs, expr.rhs]
            elif isinstance(expr, L.NaryOp):
                return expr.args
            elif isinstance(expr, L.Conditional):
                return [expr.condition, expr.true, expr.false]
            elif isinstance(expr, L.Call):
                return [expr.arguments]
            elif isinstance(expr, L.ArrayAccess):
                return [expr.array, *expr.indices]
            elif isinstance(expr, L.New):
                return []

            raise RuntimeError("Unsupported object in expression.")

        def dfs(expr: L.CExpr):
            """Depth-first search generator for CNodes expressions."""
            yield expr
            for child in child_nodes(expr):
                for sub_expr in dfs(child):
                    yield sub_expr

        def is_function_call_with_vec_args(expr: L.CExpr) -> bool:
            """Returns whether the given CNodes expression is a function call with vectorized arguments."""
            if isinstance(expr, L.Call):
                return any(was_vectorized(arg) for arg in expr.arguments)
            else:
                return False

        def is_vectorized_variable(expr: L.CExpr) -> bool:
            """Returns whether the given CNodes expression is a vectorized variable."""
            if isinstance(expr, L.Symbol) or isinstance(expr, L.ArrayAccess):
                return was_vectorized(expr)
            else:
                return False

        def is_vector_expression(expr: L.CExpr) -> bool:
            """Returns whether a CNodes expression is of vector type."""

            return ((not any(is_function_call_with_vec_args(sub_expr) for sub_expr in dfs(expr)))
                    and any(is_vectorized_variable(sub_expr) for sub_expr in dfs(expr)))

        @functools.singledispatch
        def vectorize(stmnt):
            """Overloaded function to vectorize any CNodes statement, default implementation."""

            ignore_list = [
                L.Comment,
                L.CExprLiteral,
                L.Pragma,
                L.Using,
                L.Break,
                L.Continue,
                L.Return
            ]

            if any(isinstance(stmnt, ignored_type) for ignored_type in ignore_list):
                return stmnt

            # For dev: check whether to add type to ignore list or to implement it below
            raise RuntimeError("Vectorization of {} CNodes statement not implemented!".format(type(stmnt)))

        # Vectorization of statements

        @vectorize.register(L.VariableDecl)
        def vectorize_variable_decl(stmnt):
            # Skip static variables
            if "static" in stmnt.typename:
                return stmnt

            # Vectorize the type of the variable
            var_decl = copy(stmnt)
            var_decl.typename = stmnt.typename.replace(base_type, vector_type)

            # Log that the scalar was vectorized
            ctx["vectorized_intermediates"].append(stmnt.symbol.name)

            # Check value, if it isn't a vector expression the assignment has to be wrapped in a loop
            if stmnt.value is not None:
                if not is_vector_expression(stmnt.value):
                    assign_op = L.Assign(var_decl.symbol, stmnt.value)
                    var_decl.value = None
                    return L.StatementList([var_decl, vectorize(assign_op)])

            return var_decl

        @vectorize.register(L.ArrayDecl)
        def vectorize_array_decl(stmnt):
            # Skip static arrays
            if "static" in stmnt.typename:
                return stmnt

            if stmnt.values is not None:
                if not isinstance(stmnt.values, int) and not isinstance(stmnt.values, float):
                    raise RuntimeError(
                        "Assignment of expressions in a ArrayDecl unsupported during vectorization")

            array_decl = copy(stmnt)
            array_decl.typename = stmnt.typename.replace(base_type, vector_type)

            # Store that the array was vectorized
            ctx["vectorized_intermediates"].append(stmnt.symbol.name)

            return array_decl

        @vectorize.register(L.AssignOp)
        def vectorize_assign(stmnt):
            # Transform all assignment operations (=, +=, *=,...)
            target, value = stmnt.lhs, stmnt.rhs

            if was_vectorized(target) and not is_vector_expression(value):
                assign_op = copy(stmnt)
                assign_op.lhs = vectorize(target)
                assign_op.rhs = vectorize(value)
                return simd_loop(assign_op)
            else:
                return stmnt

        @vectorize.register(L.ForRange)
        def vectorize_for_range(stmnt):
            # Vectorize for-loop by vectorizing its body
            for_range = copy(stmnt)
            for_range.body = vectorize(stmnt.body)
            return for_range

        @vectorize.register(L.Statement)
        def vectorize_statement(stmnt):
            # Statement is used to wrap assignment expressions into a statement
            stmnt = copy(stmnt)

            # Vectorize the contained expression
            vectorized_expr = vectorize(stmnt.expr)

            # The vectorized version might already be a statement
            try:
                vectorized_stmnt = L.Statement(vectorized_expr)
                return vectorized_stmnt
            except RuntimeError:
                return vectorized_expr

        @vectorize.register(L.StatementList)
        def vectorize_statement_list(stmnts):
            return L.StatementList([vectorize(stmnt) for stmnt in stmnts.statements])

        # Vectorization of expressions, only applied inside of for-loops generated by a statement above

        @vectorize.register(L.Symbol)
        def vectorize_symbol(expr):
            if was_vectorized(expr):
                return L.ArrayAccess(expr, (i_simd,))
            else:
                return expr

        @vectorize.register(L.ArrayAccess)
        def vectorize_array_access(expr):
            if was_vectorized(expr):
                array_access = copy(expr)
                array_access.indices = expr.indices + (i_simd,)
                return array_access
            else:
                return expr

        @vectorize.register(L.BinOp)
        def vectorize_binop(expr):
            # Transform all binary operators (+, *, -, %,...)
            bin_op = copy(expr)
            bin_op.lhs = vectorize(expr.lhs)
            bin_op.rhs = vectorize(expr.rhs)
            return bin_op

        @vectorize.register(L.NaryOp)
        def vectorize_sum(expr):
            # Transform all n-ary operators (sum, product,...)
            nary_op = copy(expr)
            nary_op.args = [vectorize(arg) for arg in expr.args]
            return nary_op

        @vectorize.register(L.Call)
        def vectorize_call(expr):
            # Transform functions calls, assuming that the call has no side effects
            if expr.arguments is not None:
                call = copy(expr)
                call.arguments = [vectorize(arg) for arg in expr.arguments]
                return call
            else:
                return expr

        vectorized = [vectorize(stmnt) for stmnt in statements]
        return vectorized

    def vectorize_with_loops(self, statements, vec_length=4, alignment=32):
        """
        Converts a list of tabulate_tensor CNodes statements into a cross-element vectorized version.
        :param statements: The list of statements that should be transformed
        :param vec_length: The number of elements to perform cross element vectorization over
        :param alignment: Alignment in bytes used for arrays whose rank is increased for vectorization
        :return: A list of the transformed statements
        """

        L = self.backend.language

        enable_cross_element_fuse = self.ir["params"]["enable_cross_element_fuse"]
        enable_cross_element_array_conv = self.ir["params"]["enable_cross_element_array_conv"]

        ctx = {
            "vectorized_inputs": ["A", "w", "coordinate_dofs"],
            "vectorized_intermediates": [],
            "reduce_to_scalars": [],
            "reduced_typenames": dict(),
            "reduced_scalar_names": set()
        }

        if enable_cross_element_array_conv:
            ctx["reduce_to_scalars"].append("sp")

        # Symbol used as index in for loops over the elements
        i_simd = L.Symbol("i_elem")

        def should_reduce(expr: L.Symbol) -> bool:
            """Returns whether the symbol belongs to an array that should be reduced to a scalar."""

            return expr.name in ctx["reduce_to_scalars"]

        def was_reduced(expr: L.ArrayAccess) -> bool:
            """Returns whether the array used in this array access was reduced to a scalar."""

            scalar_name = expr.array.name + str(expr.indices[0])
            return should_reduce(expr.array) and scalar_name in ctx["reduced_scalar_names"]

        def was_expanded(expr) -> bool:
            """Returns whether an array's/scalar's rank was increased for vectorization."""

            if isinstance(expr, L.ArrayAccess):
                return expr.array.name in ctx["vectorized_intermediates"] or expr.array.name in ctx[
                    "vectorized_inputs"]
            if isinstance(expr, L.Symbol):
                return expr.name in ctx["vectorized_intermediates"]

            raise RuntimeError("Unsupported expression")

        def simd_loop(body):
            """Returns a ForRange that is used to perform a cross element loop."""
            return L.ForRange(i_simd, 0, vec_length, body)

        @functools.singledispatch
        def vectorize(stmnt):
            """Overloaded function to vectorize any CNodes statement, default implementation"""

            ignore_list = [
                L.Comment,
                L.CExprLiteral,
                L.Pragma,
                L.Using,
                L.Break,
                L.Continue,
                L.Return
            ]

            if any(isinstance(stmnt, ignored_type) for ignored_type in ignore_list):
                return stmnt

            # For dev: check whether to add type to ignore list or to implement it below
            raise RuntimeError("Vectorization of {} CNodes statement not implemented!".format(type(stmnt)))

        @vectorize.register(L.ForRange)
        def vectorize_for_range(stmnt):
            # Vectorize for loop by vectorizing its body
            vec_stmnt = copy(stmnt)
            vec_stmnt.body = vectorize(stmnt.body)
            return vec_stmnt

        @vectorize.register(L.VariableDecl)
        def vectorize_variable_decl(stmnt):
            # Skip static variables
            if "static" in stmnt.typename:
                return stmnt

            # Remove constness
            typename = stmnt.typename
            if "const" in typename:
                typename = typename.replace("const", "").strip()

            array_decl = L.ArrayDecl(typename=typename,
                                     symbol=stmnt.symbol,
                                     sizes=vec_length,
                                     alignas=alignment)

            # Log that the scalar was transformed to an array
            ctx["vectorized_intermediates"].append(stmnt.symbol.name)

            if stmnt.value is not None:
                # If the scalar had a value, it has to be assigned in a loop and vectorized recursively
                value = vectorize(stmnt.value)
                assignment = L.Assign(stmnt.symbol[i_simd], value)

                return L.StatementList([array_decl, simd_loop(assignment)])
            else:
                return array_decl

        @vectorize.register(L.ArrayDecl)
        def vectorize_array_decl(stmnt):
            # Skip static arrays
            if "static" in stmnt.typename:
                return stmnt

            if stmnt.values is not None and stmnt.values != 0:
                raise RuntimeError("Values in vectorization for ArrayDecl unsupported")

            if should_reduce(stmnt.symbol):
                ctx["reduced_typenames"][stmnt.symbol.name] = stmnt.typename
                return L.NoOp()

            array_decl = copy(stmnt)
            # Increase rank
            array_decl.sizes = stmnt.sizes + (vec_length,)

            # Log that the rank was increased
            ctx["vectorized_intermediates"].append(stmnt.symbol.name)

            return array_decl

        @vectorize.register(L.ArrayAccess)
        def vectorize_array_access(stmnt):
            if was_reduced(stmnt):
                return vectorize(L.Symbol(stmnt.array.name + str(stmnt.indices[0])))

            # Check whether rank was increased
            if was_expanded(stmnt):
                array_access = copy(stmnt)

                # Update array indexing
                if stmnt.array.name in ctx["vectorized_intermediates"]:
                    # For manually expanded arrays, simply append array index
                    array_access.indices = stmnt.indices + (i_simd,)
                elif stmnt.array.name in ctx["vectorized_inputs"]:
                    # For in/out arrays, we have strided access instead
                    array_access.indices = stmnt.indices[0:-1] + (
                        i_simd + L.LiteralInt(vec_length) * stmnt.indices[-1],)

                return array_access

            else:
                # Rank is unchanged (e.g. static array)
                return stmnt

        @vectorize.register(L.Symbol)
        def vectorize_symbol(stmnt):
            if was_expanded(stmnt):
                # Scalar which was transformed to an array:
                return L.ArrayAccess(stmnt, (i_simd,))
            else:
                return stmnt

        @vectorize.register(L.AssignOp)
        def vectorize_assign(stmnt):
            # Transform all assignment operations (=, +=, *=,...)
            target, value = stmnt.lhs, stmnt.rhs

            if isinstance(target, L.ArrayAccess) and should_reduce(target.array):
                assert len(target.indices) == 1
                new_target_name = target.array.name + str(target.indices[0])

                if was_reduced(target):
                    target = L.Symbol(new_target_name)
                else:
                    assert isinstance(stmnt, L.Assign)
                    var_decl = L.VariableDecl(ctx["reduced_typenames"][target.array.name],
                                              L.Symbol(new_target_name),
                                              value=value)
                    ctx["reduced_scalar_names"].add(new_target_name)
                    return vectorize(var_decl)

            if was_expanded(target):
                assignment = type(stmnt)(vectorize(target), vectorize(value))
                return simd_loop(assignment)
            else:
                return stmnt

        @vectorize.register(L.BinOp)
        def vectorize_binop(stmnt):
            # Transform all binary operators (+, *, -, %,...)
            return type(stmnt)(vectorize(stmnt.lhs), vectorize(stmnt.rhs))

        @vectorize.register(L.NaryOp)
        def vectorize_sum(stmnt):
            # Transform all n-ary operators (sum, product,...)
            return type(stmnt)([vectorize(arg) for arg in stmnt.args])

        @vectorize.register(L.Call)
        def vectorize_call(stmnt):
            # Transform functions calls, assuming that the call has no side effects
            if stmnt.arguments is not None:
                return L.Call(stmnt.function, [vectorize(arg) for arg in stmnt.arguments])
            else:
                return stmnt

        @vectorize.register(L.Statement)
        def vectorize_statement(stmnt):
            # Vectorize the contained expression
            vectorized_expr = vectorize(stmnt.expr)

            # The vectorized version might already be a statement
            try:
                vectorized_stmnt = L.Statement(vectorized_expr)
                return vectorized_stmnt
            except RuntimeError:
                return vectorized_expr

        @vectorize.register(L.StatementList)
        def vectorize_statement_list(stmnts):
            return L.StatementList([vectorize(stmnt) for stmnt in stmnts.statements])

        def optimize(stmnts):
            """Joins consecutive cross-element expanded assignment loops to a single loop."""

            if len(stmnts) < 2:
                return stmnts

            def get_statements_as_list(stmnt):
                """For a StatementList, return statements list, otherwise return stmnt as list."""

                if isinstance(stmnt, L.StatementList):
                    stmnts = stmnt.statements
                else:
                    stmnts = [stmnt]

                return stmnts

            def join_for_ranges(range1, range2):
                """Joins two ForRange objects."""

                assert isinstance(range1, L.ForRange)
                assert isinstance(range2, L.ForRange)
                attributes = ("index", "begin", "end", "index_type")
                assert all(getattr(range1, name) == getattr(range2, name) for name in attributes)

                range3 = copy(range1)

                stmnts1 = get_statements_as_list(range1.body)
                stmnts2 = get_statements_as_list(range2.body)

                range3.body = L.StatementList(stmnts1 + stmnts2)
                return range3

            def is_expanded_assignment(stmnt):
                """Return whether the specified statement is a cross-element expanded assignment loop."""
                return (isinstance(stmnt, L.ForRange)
                        and (stmnt.index == i_simd)
                        and (stmnt.begin.value == 0)
                        and (stmnt.end.value == vec_length))

            def is_expanded_var_decl(stmnt):
                """
                Returns whether the specified statement is a StatementList,
                that declares and assigns a cross-element expanded scalar.
                """
                if isinstance(stmnt, L.StatementList) and len(stmnt.statements) == 2:
                    decl = stmnt.statements[0]
                    assign = stmnt.statements[1]

                    if isinstance(decl, L.ArrayDecl):
                        return decl.sizes[-1] == vec_length and is_expanded_assignment(assign)

                return False

            # "Enum" for the types of statements that are recognized by optimizer
            UNINTERSTING_TYPE = 0
            EXPANDED_ASSIGNMENT = 1
            EXPANDED_VAR_DECL = 2

            def stmnt_type(stmnt):
                if is_expanded_assignment(stmnt):
                    return EXPANDED_ASSIGNMENT
                if is_expanded_var_decl(stmnt):
                    return EXPANDED_VAR_DECL

                return UNINTERSTING_TYPE

            optimized = []

            prev_stmnt = stmnts[0]
            prev_type = stmnt_type(prev_stmnt)

            # Loop over all statements in list
            for stmnt in stmnts[1:]:
                curr_type = stmnt_type(stmnt)

                # We can only join if consective statements are of same type
                if prev_type == curr_type:
                    # Expanded assignment loops can be fused
                    if curr_type == EXPANDED_ASSIGNMENT:
                        prev_stmnt = join_for_ranges(prev_stmnt, stmnt)

                        continue

                    # Expanded variable declarations can be rearranged and fused
                    elif curr_type == EXPANDED_VAR_DECL:
                        decls = [*prev_stmnt.statements[:-1], stmnt.statements[0]]
                        assign = join_for_ranges(prev_stmnt.statements[-1], stmnt.statements[1])

                        prev_stmnt = L.StatementList(decls + [assign])

                        continue

                # Otherwise, replace the previous stored statement used for comparison
                optimized.append(prev_stmnt)
                prev_stmnt = stmnt
                prev_type = curr_type

            # Append last statement
            optimized.append(prev_stmnt)
            return optimized

        vectorized = [vectorize(stmnt) for stmnt in statements]

        if enable_cross_element_fuse:
            vectorized = optimize(vectorized)

        return vectorized

    def generate_quadrature_tables(self):
        """Generate static tables of quadrature points and weights."""
        L = self.backend.language

        parts = []

        # No quadrature tables for custom (given argument)
        # or point (evaluation in single vertex)
        skip = custom_integral_types + point_integral_types
        if self.ir["integral_type"] in skip:
            return parts

        alignas = self.ir["params"]["alignas"]

        # Loop over quadrature rules
        for num_points in self.ir["all_num_points"]:
            varying_ir = self.ir["varying_irs"][num_points]

            points, weights = self.ir["quadrature_rules"][num_points]
            assert num_points == len(weights)
            assert num_points == points.shape[0]

            # Generate quadrature weights array
            if varying_ir["need_weights"]:
                wsym = self.backend.symbols.weights_table(num_points)
                parts += [
                    L.ArrayDecl("static const double", wsym, num_points, weights, alignas=alignas)
                ]

            # Generate quadrature points array
            N = product(points.shape)
            if varying_ir["need_points"] and N:
                # Flatten array: (TODO: avoid flattening here, it makes padding harder)
                flattened_points = points.reshape(N)
                psym = self.backend.symbols.points_table(num_points)
                parts += [
                    L.ArrayDecl("static const double", psym, N, flattened_points, alignas=alignas)
                ]

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
        inline_tables = self.ir["integral_type"] == "cell"

        alignas = self.ir["params"]["alignas"]
        padlen = self.ir["params"]["padlen"]

        if self.ir["integral_type"] in custom_integral_types:
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
                "static const double", name, table.shape, table, alignas=alignas, padlen=p)
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

        assert self.ir["integral_type"] in custom_integral_types

        num_points = self.ir["fake_num_points"]
        chunk_size = self.ir["params"]["chunk_size"]

        gdim = self.ir["geometric_dimension"]

        alignas = self.ir["params"]["alignas"]
        # padlen = self.ir["params"]["padlen"]

        tables = self.ir["unique_tables"]
        table_types = self.ir["unique_table_types"]
        # table_origins = self.ir["unique_table_origins"]  # FIXME

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
            varying_ir = self.ir["varying_irs"][num_points]

            # Copy quadrature weights for this chunk
            if varying_ir["need_weights"]:
                cwsym = self.backend.symbols.custom_quadrature_weights()
                wsym = self.backend.symbols.custom_weights_table()
                rule_parts += [
                    L.ArrayDecl("double", wsym, chunk_size, 0, alignas=alignas),
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
                    L.ArrayDecl("double", psym, chunk_size * gdim, 0, alignas=alignas),
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
                    "double", name, (1, chunk_size, table.shape[2]), 0,
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

        num_points = None
        expr_ir = self.ir["piecewise_ir"]

        name = "sp"
        arraysymbol = L.Symbol(name)
        parts = self.generate_partition(arraysymbol, expr_ir["V"], expr_ir["V_active"],
                                        expr_ir["V_mts"], expr_ir["mt_tabledata"], num_points)
        parts = L.commented_code_list(parts, "Unstructured piecewise computations")
        return parts

    def generate_unstructured_varying_partition(self, num_points):
        L = self.backend.language

        expr_ir = self.ir["varying_irs"][num_points]

        name = "sv"
        arraysymbol = L.Symbol("%s%d" % (name, num_points))
        parts = self.generate_partition(arraysymbol, expr_ir["V"], expr_ir["V_varying"],
                                        expr_ir["V_mts"], expr_ir["mt_tabledata"], num_points)
        parts = L.commented_code_list(parts, "Unstructured varying computations for num_points=%d" %
                                      (num_points, ))
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
                vaccess = self.backend.access(mt.terminal, mt, tabledata, num_points)
                vdef = self.backend.definitions(mt.terminal, mt, tabledata, num_points, vaccess)

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
                # j = variable_id[i]

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
                    if self.ir["params"]["use_symbol_array"]:
                        vaccess = symbol[j]
                        intermediates.append(L.Assign(vaccess, vexpr))
                    else:
                        vaccess = L.Symbol("%s_%d" % (symbol.name, j))
                        intermediates.append(L.VariableDecl("const double", vaccess, vexpr))
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
            if self.ir["params"]["use_symbol_array"]:
                alignas = self.ir["params"]["alignas"]
                parts += [L.ArrayDecl("double", symbol, len(intermediates), alignas=alignas)]
            parts += intermediates
        return parts

    def generate_dofblock_partition(self, num_points):
        if num_points is None:  # NB! None meaning piecewise partition, not custom integral
            block_contributions = self.ir["piecewise_ir"]["block_contributions"]
        else:
            block_contributions = self.ir["varying_irs"][num_points]["block_contributions"]

        preparts = []
        quadparts = []
        postparts = []

        blocks = [(blockmap, blockdata)
                  for blockmap, contributions in sorted(block_contributions.items())
                  for blockdata in contributions if blockdata.block_mode != "preintegrated"]

        for blockmap, blockdata in blocks:
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

    def get_entities(self, blockdata):
        L = self.backend.language

        if self.ir["integral_type"] == "interior_facet":
            # Get the facet entities
            entities = []
            for r in blockdata.restrictions:
                if r is None:
                    entities.append(0)
                else:
                    entities.append(self.backend.symbols.entity(self.ir["entitytype"], r))
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
                entity = self.backend.symbols.entity(self.ir["entitytype"], None)
            return (entity, )

    def get_arg_factors(self, blockdata, block_rank, num_points, iq, indices):
        L = self.backend.language

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

            table = self.backend.symbols.element_table(td, self.ir["entitytype"], mt.restriction)

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

        fwtempname = "fw"
        tempname = tempnames.get(blockdata.block_mode)

        alignas = self.ir["params"]["alignas"]
        padlen = self.ir["params"]["padlen"]

        block_rank = len(blockmap)
        blockdims = tuple(len(dofmap) for dofmap in blockmap)
        padded_blockdims = pad_innermost_dim(blockdims, padlen)

        ttypes = blockdata.ttypes
        if "zeros" in ttypes:
            raise FFCError("Not expecting zero arguments to be left in dofblock generation.")

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
            preparts.append(L.ArrayDecl("double", B, blockdims, 0, alignas=alignas, padlen=padlen))

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
                fw, defined = self.get_temp_symbol(fwtempname, key)
                if not defined:
                    quadparts.append(L.VariableDecl("const double", fw, fw_rhs))

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
                        L.ArrayDecl("double", P, P_dim, None, alignas=alignas, padlen=padlen))
                    P_rhs = L.float_product([fw, arg_factors[i]])
                    body = L.Assign(P[P_index], P_rhs)
                    # if ttypes[i] != "quadrature":  # FIXME: What does this mean here?
                    vectorize = self.ir["params"]["vectorize"]
                    body = L.ForRange(P_index, 0, P_dim, body=body, vectorize=vectorize)
                    quadparts.append(body)

                B_rhs = P[P_index] * arg_factors[j]

            # Add result to block inside quadloop
            body = L.AssignAdd(B[B_indices], B_rhs)  # NB! += not =
            for i in reversed(range(block_rank)):
                # Vectorize only the innermost loop
                vectorize = self.ir["params"]["vectorize"] and (i == block_rank - 1)
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
                preparts.append(L.ArrayDecl("double", P, P_dim, 0, alignas=alignas, padlen=padlen))

                # Multiply collected factors
                P_rhs = L.float_product([fw, arg_factors[not_piecewise_index]])

                # Accumulate P += weight * f * args in quadrature loop
                body = L.AssignAdd(P[P_index], P_rhs)
                body = L.ForRange(P_index, 0, pad_dim(P_dim, padlen), body=body)
                quadparts.append(body)

            # Define B = B_rhs = piecewise_argument[:] * P[:], where P[:] = sum_q weight * f * other_argument[:]
            B_rhs = arg_factors[i] * P[P_index]

            # Define rhs expression for A[blockmap[arg_indices]] += A_rhs
            A_rhs = B_rhs

        elif blockdata.block_mode in ("premultiplied", "preintegrated"):
            P_entity_indices = self.get_entities(blockdata)
            if blockdata.transposed:
                P_block_indices = (arg_indices[1], arg_indices[0])
            else:
                P_block_indices = arg_indices
            P_ii = P_entity_indices + P_block_indices

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
                    preparts += [L.VariableDecl("double", FI, 0)]
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
        # FIXME: Generalize this to unrolling all A[] += ... loops, or all loops with noncontiguous DM??
        L = self.backend.language

        block_contributions = self.ir["piecewise_ir"]["block_contributions"]

        blocks = [(blockmap, blockdata)
                  for blockmap, contributions in sorted(block_contributions.items())
                  for blockdata in contributions if blockdata.block_mode == "preintegrated"]

        # Get symbol, dimensions, and loop index symbols for A
        A_shape = self.ir["tensor_shape"]
        A_size = product(A_shape)
        A_rank = len(A_shape)

        # TODO: there's something like shape2strides(A_shape) somewhere
        A_strides = [1] * A_rank
        for i in reversed(range(0, A_rank - 1)):
            A_strides[i] = A_strides[i + 1] * A_shape[i + 1]

        A_values = [0.0] * A_size

        for blockmap, blockdata in blocks:
            # Accumulate A[blockmap[...]] += f*PI[...]

            # Get table for inlining
            tables = self.ir["unique_tables"]
            table = tables[blockdata.name]
            inline_table = self.ir["integral_type"] == "cell"

            # Get factor expression
            v = self.ir["piecewise_ir"]["V"][blockdata.factor_index]
            f = self.get_var(None, v)

            # Define rhs expression for A[blockmap[arg_indices]] += A_rhs
            # A_rhs = f * PI where PI = sum_q weight * u * v
            PI = L.Symbol(blockdata.name)
            # block_rank = len(blockmap)

            # # Override dof index with quadrature loop index for arguments with
            # # quadrature element, to index B like B[iq*num_dofs + iq]
            # arg_indices = tuple(
            #     self.backend.symbols.argument_loop_index(i) for i in range(block_rank))

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

        init_mode = self.ir["params"]["tensor_init_mode"]
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
            raise FFCError("Invalid init_mode parameter %s" % (init_mode, ))

        return parts

    def generate_expr_copyout_statements(self):
        L = self.backend.language
        parts = []

        # Not expecting any quadrature loop scopes here
        assert tuple(self.scopes.keys()) == (None, )

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

        return parts

    def generate_tensor_copyout_statements(self):
        L = self.backend.language
        parts = []

        # Get symbol, dimensions, and loop index symbols for A
        A_shape = self.ir["tensor_shape"]
        A_rank = len(A_shape)

        # TODO: there's something like shape2strides(A_shape) somewhere
        A_strides = [1] * A_rank
        for i in reversed(range(0, A_rank - 1)):
            A_strides[i] = A_strides[i + 1] * A_shape[i + 1]

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

    def generate_copyout_statements(self):
        """Generate statements copying results to output array."""
        if self.ir["integral_type"] == "expression":
            return self.generate_expr_copyout_statements()
        # elif self.ir["unroll_copyout_loops"]:
        #    return self.generate_unrolled_tensor_copyout_statements()
        else:
            return self.generate_tensor_copyout_statements()
