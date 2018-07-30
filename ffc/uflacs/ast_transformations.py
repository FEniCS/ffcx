# -*- coding: utf-8 -*-
# Copyright (C) 2018 Fabian Löschner
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Algorithms that perform transformation of the AST of the tabulate_tensor functions"""

import functools
from copy import copy


class Vectorizer(object):
    def __init__(self, ir, backend):
        self.backend = backend
        self.ir = ir

        self.vec_length = self.ir["integrals_metadata"]["cell_batch_size"]
        self.align = self.ir["params"]["alignas"]

    def vectorize(self, statements):
        # Optionally use gcc vector extensions
        if self.ir["params"]["enable_cross_element_gcc_ext"]:
            vectorized = self.__vectorize_with_gcc_exts(statements)
        else:
            vectorized = self.__vectorize_with_loops(statements)

        return vectorized

    def __vectorize_with_gcc_exts(self, statements):
        """Cross-element vectorization of `tabulate_tensor` statement list with GCC extensions.

        Converts a list of `tabulate_tensor` CNodes statements into a cross-element vectorized version using
        GCC's vector extensions.

        :param statements: The list of statements that should be transformed
        :param vec_length: The number of elements to perform cross element vectorization over
        :param alignment: Alignment in bytes used for arrays whose rank is increased for vectorization
        :return: A list of the transformed statements
        """

        L = self.backend.language

        base_type = "double"
        vec_length = self.vec_length
        vector_type = "{}{}".format(base_type, str(vec_length))

        vectorized_parameters = self.__get_cell_dependent_parameters()

        ctx = {
            "vectorized_parameters": [param_decl.symbol.name for param_decl in vectorized_parameters],
            "vectorized_intermediates": [],
        }

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
                (symbol_name in arr) for arr in [ctx["vectorized_intermediates"], ctx["vectorized_parameters"]])

        def dfs(cnode: L.CNode):
            """Depth-first search generator for CNodes expressions."""
            yield cnode
            for child in cnode.children():
                for childs_child in dfs(child):
                    yield childs_child

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

        # Run vectorization of all statements
        vectorized = [vectorize(stmnt) for stmnt in statements]

        # ---------
        # Handle casting of function parameters to vector extension types

        preamble = [L.VerbatimStatement("typedef double {} __attribute__ ((vector_size ({})));".format(vector_type,
                                                                                                       8 * vec_length))]

        def get_vectorized_name(name: str) -> str:
            if name in ctx["vectorized_parameters"]:
                return name + "_x"
            return name

        # Keep track of actually used parameters to only create cast statements for them
        used_parameters = set()

        # Modify all symbol names that refer to parameters
        for stmnt in vectorized:
            for child in dfs(stmnt):
                if isinstance(child, L.Symbol):
                    if child.name in ctx["vectorized_parameters"]:
                        used_parameters.add(child.name)
                        child.name = get_vectorized_name(child.name)

        # Generate casts to vector type of all used parameters
        for param_decl in vectorized_parameters:
            old_name = param_decl.symbol.name

            # Skip parameters that were not used
            if old_name not in used_parameters:
                continue

            new_name = get_vectorized_name(old_name)
            new_typename = param_decl.typename.replace(base_type, vector_type)

            # Append the cast statement
            cast = L.VerbatimStatement("{0} {1} = ({0}){2};".format(new_typename, new_name, old_name))
            preamble.append(cast)

        # Join vectorized code with cast statements
        vectorized = preamble + vectorized
        return vectorized

    def __vectorize_with_loops(self, statements):
        """Cross-element vectorization of `tabulate_tensor` statement list.

        Converts a list of `tabulate_tensor` CNodes statements into a cross-element vectorized version.

        :param statements: The list of statements that should be transformed
        :param vec_length: The number of elements to perform cross element vectorization over
        :param alignment: Alignment in bytes used for arrays whose rank is increased for vectorization
        :return: A list of the transformed statements
        """

        L = self.backend.language

        vec_length = self.vec_length
        align = self.align

        enable_cross_element_fuse = self.ir["params"]["enable_cross_element_fuse"]
        enable_cross_element_array_conv = self.ir["params"]["enable_cross_element_array_conv"]

        ctx = {
            "vectorized_parameters": ["A", "w", "coordinate_dofs"],
            "vectorized_intermediates": [],
            "reduce_to_scalars": [],
            "reduced_typenames": {},
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

            # FIXME: why not expr.indices[-1]?
            # Catch the case where a scalar was added as an array ot the AST
            index = expr.indices[0] if len(expr.indices) > 0 else 0
            scalar_name = expr.array.name + str(index)
            return should_reduce(expr.array) and scalar_name in ctx["reduced_scalar_names"]

        def was_expanded(expr) -> bool:
            """Returns whether an array's/scalar's rank was increased for vectorization."""

            if isinstance(expr, L.ArrayAccess):
                return expr.array.name in ctx["vectorized_intermediates"] or expr.array.name in ctx[
                    "vectorized_parameters"]
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
                                     alignas=align)

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
                elif stmnt.array.name in ctx["vectorized_parameters"]:
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
                """Checks if statement is a cross-element version of an AssignOp.

                Returns whether the specified statement is a cross-element expanded assignment loop.
                """
                return (isinstance(stmnt, L.ForRange)
                        and (stmnt.index == i_simd)
                        and (stmnt.begin.value == 0)
                        and (stmnt.end.value == vec_length))

            def is_expanded_var_decl(stmnt):
                """Checks if statement is a cross-element version of a VarDecl.

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

    def __get_cell_dependent_parameters(self):
        L = self.backend.language

        inputs = {
            "cell": [
                L.VariableDecl("double* restrict", "A"),
                L.VariableDecl("const double* const*", "w"),
                L.VariableDecl("const double* restrict", "coordinate_dofs"),
            ]
        }

        return inputs[self.ir["integral_type"]]
