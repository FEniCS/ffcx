# Copyright (C) 2025 Chris Richardson and Paul T. Kühner
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Numba implementation for output."""

from functools import singledispatchmethod

import numpy as np
from numpy import typing as npt

import ffcx.codegeneration.lnodes as L
from ffcx.codegeneration.interface import Formatter as FormatterInterface
from ffcx.codegeneration.utils import dtype_to_scalar_dtype


def build_initializer_lists(values: npt.NDArray) -> str:
    """Build list of values."""
    arr = "["
    if len(values.shape) == 1:
        return "[" + ", ".join(str(v) for v in values) + "]"
    elif len(values.shape) > 1:
        arr += ",\n".join(build_initializer_lists(v) for v in values)
    arr += "]"
    return arr


class Formatter(FormatterInterface):
    """Implementation for numba output backend."""

    scalar_type: np.dtype
    real_type: np.dtype

    def __init__(self, dtype: npt.DTypeLike) -> None:
        """Initialise."""
        self.scalar_type = np.dtype(dtype)
        self.real_type = dtype_to_scalar_dtype(dtype)

    def _dtype_to_name(self, dtype: L.DataType) -> str:
        """Convert dtype to Python name."""
        if dtype == L.DataType.SCALAR:
            return f"np.{self.scalar_type}"
        if dtype == L.DataType.REAL:
            return f"np.{self.real_type}"
        if dtype == L.DataType.INT:
            return f"np.{np.int32}"
        if dtype == L.DataType.BOOL:
            return f"np.{np.bool}"
        raise ValueError(f"Invalid dtype: {dtype}")

    @singledispatchmethod
    def __call__(self, obj: L.LNode) -> str:
        """Format an L Node."""
        raise NotImplementedError(f"Can not format object to type {type(obj)}")

    @__call__.register
    def _(self, section: L.Section) -> str:
        """Format a section."""
        # add new line before section
        comments = self._format_comment_str("------------------------")
        comments += self._format_comment_str(f"Section: {section.name}")
        comments += self._format_comment_str(f"Inputs: {', '.join(w.name for w in section.input)}")
        comments += self._format_comment_str(
            f"Outputs: {', '.join(w.name for w in section.output)}"
        )
        declarations = "".join(self(s) for s in section.declarations)

        body = ""
        if len(section.statements) > 0:
            body = "".join(self(s) for s in section.statements)

        body += self._format_comment_str("------------------------")
        return comments + declarations + body

    @__call__.register
    def _(self, slist: L.StatementList) -> str:
        """Format a list of statements."""
        output = ""
        for s in slist.statements:
            output += self(s)
        return output

    def _format_comment_str(self, comment: str) -> str:
        """Format str to comment string."""
        return f"# {comment} \n"

    @__call__.register
    def _(self, c: L.Comment) -> str:
        """Format a comment."""
        return self._format_comment_str(c.comment)

    @__call__.register
    def _(self, arr: L.ArrayDecl) -> str:
        """Format an array declaration."""
        dtype = arr.symbol.dtype
        typename = self._dtype_to_name(dtype)

        symbol = self(arr.symbol)
        if arr.values is None:
            return f"{symbol} = np.empty({arr.sizes}, dtype={typename})\n"
        elif arr.values.size == 1:
            return f"{symbol} = np.full({arr.sizes}, {arr.values[0]}, dtype={typename})\n"
        av = build_initializer_lists(arr.values)
        av = f"np.array({av}, dtype={typename})"
        return f"{symbol} = {av}\n"

    @__call__.register
    def _(self, arr: L.ArrayAccess) -> str:
        """Format array access."""
        array = self(arr.array)
        idx = ", ".join(self(ix) for ix in arr.indices)
        return f"{array}[{idx}]"

    @__call__.register
    def _(self, index: L.MultiIndex) -> str:
        """Format a multi-index."""
        return self(index.global_index)

    @__call__.register
    def _(self, v: L.VariableDecl) -> str:
        """Format a variable declaration."""
        sym = self(v.symbol)
        val = self(v.value)
        return f"{sym} = {val}\n"

    @__call__.register
    def _(self, oper: L.NaryOp) -> str:
        """Format a n argument operation."""
        # Format children
        args = [self(arg) for arg in oper.args]

        # Apply parentheses
        for i in range(len(args)):
            if oper.args[i].precedence >= oper.precedence:
                args[i] = f"({args[i]})"

        # Return combined string
        return f" {oper.op} ".join(args)

    @__call__.register
    def _(self, oper: L.BinOp) -> str:
        """Format a binary operation."""
        # Format children
        lhs = self(oper.lhs)
        rhs = self(oper.rhs)

        # Apply parentheses
        if oper.lhs.precedence >= oper.precedence:
            lhs = f"({lhs})"
        if oper.rhs.precedence >= oper.precedence:
            rhs = f"({rhs})"

        # Return combined string
        return f"{lhs} {oper.op} {rhs}"

    @__call__.register(L.Neg)
    @__call__.register(L.Not)
    def _(self, oper: L.Not | L.Neg) -> str:
        """Format a unary operation."""
        arg = self(oper.arg)
        if oper.arg.precedence >= oper.precedence:
            return f"{oper.op}({arg})"
        return f"{oper.op}{arg}"

    @__call__.register(L.And)
    @__call__.register(L.Or)
    def _(self, oper: L.And | L.Or) -> str:
        """Format and or or operation."""
        # Format children
        lhs = self(oper.lhs)
        rhs = self(oper.rhs)

        # Apply parentheses
        if oper.lhs.precedence >= oper.precedence:
            lhs = f"({lhs})"
        if oper.rhs.precedence >= oper.precedence:
            rhs = f"({rhs})"

        opstr = {"||": "or", "&&": "and"}[oper.op]

        # Return combined string
        return f"{lhs} {opstr} {rhs}"

    @__call__.register
    def _(self, val: L.LiteralFloat) -> str:
        """Format a literal float."""
        return f"{val.value}"

    @__call__.register
    def _(self, val: L.LiteralInt) -> str:
        """Format a literal int."""
        return f"{val.value}"

    @__call__.register
    def _(self, r: L.ForRange) -> str:
        """Format a loop over a range."""
        begin = self(r.begin)
        end = self(r.end)
        index = self(r.index)
        output = f"for {index} in range({begin}, {end}):\n"
        b = self(r.body).split("\n")
        for line in b:
            output += f"    {line}\n"
        return output

    @__call__.register
    def _(self, s: L.Statement) -> str:
        """Format a statement."""
        return self(s.expr)

    @__call__.register(L.Assign)
    @__call__.register(L.AssignAdd)
    def _(self, expr: L.Assign | L.AssignAdd) -> str:
        """Format an assignment."""
        rhs = self(expr.rhs)
        lhs = self(expr.lhs)
        return f"{lhs} {expr.op} {rhs}\n"

    @__call__.register
    def _(self, s: L.Conditional) -> str:
        """Format a conditional."""
        # Format children
        c = self(s.condition)
        t = self(s.true)
        f = self(s.false)

        # Apply parentheses
        if s.condition.precedence >= s.precedence:
            c = f"({c})"
        if s.true.precedence >= s.precedence:
            t = f"({t})"
        if s.false.precedence >= s.precedence:
            f = f"({f})"

        # Return combined string
        return f"({t} if {c} else {f})"

    @__call__.register
    def _(self, s: L.Symbol) -> str:
        """Format a symbol."""
        return f"{s.name}"

    @__call__.register
    def _(self, f: L.MathFunction) -> str:
        """Format a math function."""
        function_map = {
            "ln": "log",
            "acos": "arccos",
            "asin": "arcsin",
            "atan": "arctan",
            "atan2": "arctan2",
            "acosh": "arccosh",
            "asinh": "arcsinh",
            "atanh": "arctanh",
        }
        function = function_map.get(f.function, f.function)
        args = [self(arg) for arg in f.args]
        if "bessel_y" in function:
            return "scipy.special.yn"
        if "bessel_j" in function:
            return "scipy.special.jn"
        if function == "erf":
            return f"math.erf({args[0]})"
        argstr = ", ".join(args)
        return f"np.{function}({argstr})"
