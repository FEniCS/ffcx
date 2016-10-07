# -*- coding: utf-8 -*-
"This file implements a base class to represent an expression."

# Copyright (C) 2009-2010 Kristian B. Oelgaard
#
# This file is part of FFC.
#
# FFC is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# FFC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with FFC. If not, see <http://www.gnu.org/licenses/>.
#
# First added:  2009-08-08
# Last changed: 2010-01-21

# FFC quadrature modules.
from .symbolics import create_float


class Expr(object):
    __slots__ = ("val", "t", "_prec", "_repr", "_hash")

    def __init__(self):
        """An Expr object contains:

        val     - float, holds value of object.
        t       - Type (int), one of CONST, GEO, IP, BASIS.
        _prec   - int, precedence which is used for comparison and comparing classes.
        _repr   - str, string value of __repr__(), we only compute it once.
        _hash   - int, hash value of __hash__(), we only compute it once.

        The constructor is empty, so initialisation of variables are left to
        child classes."""
        pass

    # Representation of the expression.
    def __repr__(self):
        "Representation of the expression for comparison and debugging."
        return self._repr

    # Hash.
    def __hash__(self):
        "Hash (for lookup in {})."
        return self._hash

    # Comparison.
    def __eq__(self, other):
        "==, True if representations are equal."
        if isinstance(other, Expr):
            return self._repr == other._repr
        return False

    def __ne__(self, other):
        "!=, True if representations are not equal."
        if isinstance(other, Expr):
            return self._repr != other._repr
        return True

    def __lt__(self, other):
        """<, compare precedence and _repr if two objects have the same
precedence."""
        if not isinstance(other, Expr):
            return False
        if self._prec < other._prec:
            return True
        elif self._prec == other._prec:
            return self._repr < other._repr
        return False

    def __gt__(self, other):
        ">, opposite of __lt__."
        if not isinstance(other, Expr):
            return True
        if self._prec > other._prec:
            return True
        elif self._prec == other._prec:
            return self._repr > other._repr
        return False

    # Public functions (for FloatValue, other classes should overload
    # as needed)
    def expand(self):
        """Expand the expression.
        (FloatValue and Symbol are expanded by construction)."""
        # Nothing to be done.
        return self

    def get_unique_vars(self, var_type):
        "Get unique variables (Symbols) as a set."
        # A float is not a variable.
        return set()

    def get_var_occurrences(self):
        """Determine the number of times all variables occurs in the
        expression.  Returns a dictionary of variables and the number
        of times they occur.  Works for FloatValue and Symbol.

        """
        # There is only one float value (if it is not -1 or 1).
        if self.val == 1.0 or self.val == -1.0:
            return {}
        return {self: 1}

    def ops(self):
        """Return number of operations to compute the expression.
        This is always zero for a FloatValue."""
        # Just return 0.
        # NOTE: This means that minus in e.g., -2  and -2*x is not counted.
        return 0

    def reduce_ops(self):
        """Reduce number of operations to evaluate the expression.  There is
        nothing to be done for FloatValue and Symbol.

        """
        # Nothing to be done.
        return self

    def reduce_var(self, var):
        """Reduce the expression by another variable by using division.  This
        works for FloatValue, Symbol and Product.

        """
        return self/var

    def reduce_vartype(self, var_type):
        """Reduce expression with given var_type. It returns a tuple (found,
        remain), where 'found' is an expression that only has
        variables of type == var_type. If no variables are found,
        found=(). The 'remain' part contains the leftover after
        division by 'found' such that: self = found*remain.  Works for
        FloatValue and Symbol.

        """
        if self.t == var_type:
            return [(self, create_float(1))]
        return [((), self)]
