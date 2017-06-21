# -*- coding: utf-8 -*-
"This file implements a class to represent a product."

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
# First added:  2009-07-12
# Last changed: 2010-03-11

from functools import reduce

# FFC modules
from ffc.log import error
from ffc.quadrature.cpp import format

# FFC quadrature modules
from .symbolics import create_float
from .symbolics import create_product
from .symbolics import create_sum
from .symbolics import create_fraction
from .expr import Expr

# FFC quadrature modules
from .floatvalue import FloatValue


class Product(Expr):
    __slots__ = ("vrs", "_expanded")

    def __init__(self, variables):
        """Initialise a Product object, it derives from Expr and contains
        the additional variables:

        vrs       - a list of variables
        _expanded - object, an expanded object of self, e.g.,
                    self = x*(2+y) -> self._expanded = (2*x + x*y) (a sum), or
                    self = 2*x -> self._expanded = 2*x (self).
        NOTE: self._prec = 2."""

        # Initialise value, list of variables, class.
        self.val = 1.0
        self.vrs = []
        self._prec = 2

        # Initially set _expanded to True.
        self._expanded = True

        # Process variables if we have any.
        if variables:
            # Remove nested Products and test for expansion.
            float_val = 1.0
            for var in variables:
                # If any value is zero the entire product is zero.
                if var.val == 0.0:
                    self.val = 0.0
                    self.vrs = [create_float(0.0)]
                    float_val = 0.0
                    break

                # Collect floats into one variable
                if var._prec == 0:  # float
                    float_val *= var.val
                    continue
                # Take care of product such that we don't create
                # nested products.
                elif var._prec == 2:  # prod
                    # If expanded product is a float, just add it.
                    if var._expanded and var._expanded._prec == 0:
                        float_val *= var._expanded.val
                    # If expanded product is symbol, this product is
                    # still expanded and add symbol.
                    elif var._expanded and var._expanded._prec == 1:
                        self.vrs.append(var._expanded)
                    # If expanded product is still a product, add the
                    # variables.
                    elif var._expanded and var._expanded._prec == 2:
                        # Add copies of the variables of other product
                        # (collect floats).
                        if var._expanded.vrs[0]._prec == 0:
                            float_val *= var._expanded.vrs[0].val
                            self.vrs += var._expanded.vrs[1:]
                            continue
                        self.vrs += var._expanded.vrs
                    # If expanded product is a sum or fraction, we
                    # must expand this product later.
                    elif var._expanded and var._expanded._prec in (3, 4):
                        self._expanded = False
                        self.vrs.append(var._expanded)
                    # Else the product is not expanded, and we must
                    # expand this one later
                    else:
                        self._expanded = False
                        # Add copies of the variables of other product
                        # (collect floats).
                        if var.vrs[0]._prec == 0:
                            float_val *= var.vrs[0].val
                            self.vrs += var.vrs[1:]
                            continue
                        self.vrs += var.vrs
                    continue
                # If we have sums or fractions in the variables the
                # product is not expanded.
                elif var._prec in (3, 4):  # sum or frac
                    self._expanded = False

                # Just add any variable at this point to list of new
                # vars.
                self.vrs.append(var)

            # If value is 1 there is no need to include it, unless it
            # is the only parameter left i.e., 2*0.5 = 1.
            if float_val and float_val != 1.0:
                self.val = float_val
                self.vrs.append(create_float(float_val))
            # If we no longer have any variables add the float.
            elif not self.vrs:
                self.val = float_val
                self.vrs = [create_float(float_val)]
            # If 1.0 is the only value left, add it.
            elif abs(float_val - 1.0) < format["epsilon"] and not self.vrs:
                self.val = 1.0
                self.vrs = [create_float(1)]

        # If we don't have any variables the product is zero.
        else:
            self.val = 0.0
            self.vrs = [create_float(0)]

        # The type is equal to the lowest variable type.
        self.t = min([v.t for v in self.vrs])

        # Sort the variables such that comparisons work.
        self.vrs.sort()

        # Compute the representation now, such that we can use it
        # directly in the __eq__ and __ne__ methods (improves
        # performance a bit, but only when objects are cached).
        self._repr = "Product([%s])" % ", ".join([v._repr for v in self.vrs])

        # Use repr as hash value.
        self._hash = hash(self._repr)

        # Store self as expanded value, if we did not encounter any
        # sums or fractions.
        if self._expanded:
            self._expanded = self

    # Print functions.
    def __str__(self):
        "Simple string representation which will appear in the generated code."
        # If we have more than one variable and the first float is -1
        # exlude the 1.
        if len(self.vrs) > 1 and self.vrs[0]._prec == 0 and self.vrs[0].val == -1.0:
            # Join string representation of members by multiplication
            return format["sub"](["", format["mul"]([str(v) for v in self.vrs[1:]])])
        return format["mul"]([str(v) for v in self.vrs])

    # Binary operators.
    def __add__(self, other):
        "Addition by other objects."
        # NOTE: Assuming expanded variables.
        # If two products are equal, add their float values.
        if other._prec == 2 and self.get_vrs() == other.get_vrs():
            # Return expanded product, to get rid of 3*x + -2*x -> x, not 1*x.
            return create_product([create_float(self.val + other.val)] + list(self.get_vrs())).expand()

        elif other._prec == 1:  # sym
            if self.get_vrs() == (other,):
                # Return expanded product, to get rid of -x + x -> 0,
                # not product(0).
                return create_product([create_float(self.val + 1.0),
                                       other]).expand()
        # Return sum
        return create_sum([self, other])

    def __sub__(self, other):
        "Subtract other objects."
        if other._prec == 2 and self.get_vrs() == other.get_vrs():
            # Return expanded product, to get rid of 3*x + -2*x -> x,
            # not 1*x.
            return create_product([create_float(self.val - other.val)] + list(self.get_vrs())).expand()
        elif other._prec == 1:  # sym
            if self.get_vrs() == (other,):
                # Return expanded product, to get rid of -x + x -> 0,
                # not product(0).
                return create_product([create_float(self.val - 1.0),
                                       other]).expand()
        # Return sum
        return create_sum([self, create_product([FloatValue(-1), other])])

    def __mul__(self, other):
        "Multiplication by other objects."
        # If product will be zero.
        if self.val == 0.0 or other.val == 0.0:
            return create_float(0)

        # If other is a Sum or Fraction let them handle it.
        if other._prec in (3, 4):  # sum or frac
            return other.__mul__(self)

        # NOTE: We expect expanded sub-expressions with no nested
        # operators.  Create new product adding float or symbol.
        if other._prec in (0, 1):  # float or sym
            return create_product(self.vrs + [other])
        # Create new product adding all variables from other Product.
        return create_product(self.vrs + other.vrs)

    def __truediv__(self, other):
        "Division by other objects."
        # If division is illegal (this should definitely not happen).
        if other.val == 0.0:
            error("Division by zero.")

        # If fraction will be zero.
        if self.val == 0.0:
            return self.vrs[0]

        # If other is a Sum we can only return a fraction.
        # NOTE: Expect that other is expanded i.e., x + x -> 2*x which
        # can be handled
        # TODO: Fix x / (x + x*y) -> 1 / (1 + y).
        # Or should this be handled when reducing a fraction?
        if other._prec == 3:  # sum
            return create_fraction(self, other)

        # Handle division by FloatValue, Symbol, Product and Fraction.
        # NOTE: assuming that we get expanded variables.

        # Copy numerator, and create list for denominator.
        num = self.vrs[:]
        denom = []
        # Add floatvalue, symbol and products to the list of denominators.
        if other._prec in (0, 1):  # float or sym
            denom = [other]
        elif other._prec == 2:  # prod
            # Get copy.
            denom = other.vrs[:]
        # fraction.
        else:
            error("Did not expected to divide by fraction.")

        # Loop entries in denominator and remove from numerator (and
        # denominator).
        for d in denom[:]:
            # Add the inverse of a float to the numerator and continue.
            if d._prec == 0:  # float
                num.append(create_float(1.0 / d.val))
                denom.remove(d)
                continue
            if d in num:
                num.remove(d)
                denom.remove(d)

        # Create appropriate return value depending on remaining data.
        if len(num) > 1:
            # TODO: Make this more efficient?
            # Create product and expand to reduce
            # Product([5, 0.2]) == Product([1]) -> Float(1).
            num = create_product(num).expand()
        elif num:
            num = num[0]
        # If all variables in the numerator has been eliminated we
        # need to add '1'.
        else:
            num = create_float(1)

        if len(denom) > 1:
            return create_fraction(num, create_product(denom))
        elif denom:
            return create_fraction(num, denom[0])
        # If we no longer have a denominater, just return the
        # numerator.
        return num

    __div__ = __truediv__

    # Public functions.
    def expand(self):
        "Expand all members of the product."
        # If we just have one variable, compute the expansion of it
        # (it is not a Product, so it should be safe). We need this to
        # get rid of Product([Symbol]) type expressions.
        if len(self.vrs) == 1:
            self._expanded = self.vrs[0].expand()
            return self._expanded

        # If product is already expanded, simply return the expansion.
        if self._expanded:
            return self._expanded

        # Sort variables such that we don't call the '*' operator more
        # than we have to.
        float_syms = []
        sum_fracs = []
        for v in self.vrs:
            if v._prec in (0, 1):  # float or sym
                float_syms.append(v)
                continue
            exp = v.expand()

            # If the expanded expression is a float, sym or product,
            # we can add the variables.
            if exp._prec in (0, 1):  # float or sym
                float_syms.append(exp)
            elif exp._prec == 2:  # prod
                float_syms += exp.vrs
            else:
                sum_fracs.append(exp)
        # If we have floats or symbols add the symbols to the rest as a single
        # product (for speed).
        if len(float_syms) > 1:
            sum_fracs.append(create_product(float_syms))
        elif float_syms:
            sum_fracs.append(float_syms[0])

        # Use __mult__ to reduce list to one single variable.
        # TODO: Can this be done more efficiently without creating all
        # the intermediate variables?
        self._expanded = reduce(lambda x, y: x * y, sum_fracs)
        return self._expanded

    def get_unique_vars(self, var_type):
        "Get unique variables (Symbols) as a set."
        # Loop all members and update the set.
        var = set()
        for v in self.vrs:
            var.update(v.get_unique_vars(var_type))
        return var

    def get_var_occurrences(self):
        """Determine the number of times all variables occurs in the
        expression.  Returns a dictionary of variables and the number
        of times they occur.

        """
        # TODO: The product should be expanded at this stage, should
        # we check this?
        # Create dictionary and count number of occurrences of each
        # variable.
        d = {}
        for v in self.vrs:
            if v in d:
                d[v] += 1
                continue
            d[v] = 1
        return d

    def get_vrs(self):
        "Return all 'real' variables."
        # A product should only have one float value after
        # initialisation.
        # TODO: Use this knowledge directly in other classes?
        if self.vrs[0]._prec == 0:  # float
            return tuple(self.vrs[1:])
        return tuple(self.vrs)

    def ops(self):
        "Get the number of operations to compute product."
        # It takes n-1 operations ('*') for a product of n members.
        op = len(self.vrs) - 1

        # Loop members and add their count.
        for v in self.vrs:
            op += v.ops()

        # Subtract 1, if the first member is -1 i.e., -1*x*y -> x*y is
        # only 1 op.
        if self.vrs[0]._prec == 0 and self.vrs[0].val == -1.0:
            op -= 1
        return op

    def reduce_ops(self):
        "Reduce the number of operations to evaluate the product."
        # It's not possible to reduce a product if it is already
        # expanded and it should be at this stage.
        # TODO: Is it safe to return self.expand().reduce_ops() if
        # product is not expanded? And do we want to?

        # TODO: This should crash if it goes wrong
        return self._expanded

    def reduce_vartype(self, var_type):
        """Reduce expression with given var_type. It returns a tuple (found,
        remain), where 'found' is an expression that only has
        variables of type == var_type. If no variables are found,
        found=(). The 'remain' part contains the leftover after
        division by 'found' such that: self = found*remain.

        """
        # Sort variables according to type.
        found = []
        remains = []
        for v in self.vrs:
            if v.t == var_type:
                found.append(v)
                continue
            remains.append(v)

        # Create appropriate object for found.
        if len(found) > 1:
            found = create_product(found)
        elif found:
            found = found.pop()
        # We did not find any variables.
        else:
            return [((), self)]

        # Create appropriate object for remains.
        if len(remains) > 1:
            remains = create_product(remains)
        elif remains:
            remains = remains.pop()
        # We don't have anything left.
        else:
            return [(self, create_float(1))]

        # Return whatever we found.
        return [(found, remains)]
