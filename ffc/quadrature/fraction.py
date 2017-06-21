# -*- coding: utf-8 -*-
"This file implements a class to represent a fraction."

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
# Last changed: 2010-02-09

# FFC modules.
from ffc.log import error
from ffc.quadrature.cpp import format

# FFC quadrature modules.
from .symbolics import create_float
from .symbolics import create_product
from .symbolics import create_sum
from .symbolics import create_fraction
from .expr import Expr

# FFC quadrature modules.
from .floatvalue import FloatValue


class Fraction(Expr):
    __slots__ = ("num", "denom", "_expanded", "_reduced")

    def __init__(self, numerator, denominator):
        """Initialise a Fraction object, it derives from Expr and contains
        the additional variables:

        num       - expr, the numerator.
        denom     - expr, the denominator.
        _expanded - object, an expanded object of self, e.g.,
                    self = 'x*y/x'-> self._expanded = y (a symbol).
        _reduced  - object, a reduced object of self, e.g.,
                    self = '(2*x + x*y)/z'-> self._reduced = x*(2 + y)/z (a fraction).
        NOTE: self._prec = 4."""

        # Check for illegal division.
        if denominator.val == 0.0:
            error("Division by zero.")

        # Initialise all variables.
        self.val = numerator.val
        self.t = min([numerator.t, denominator.t])
        self.num = numerator
        self.denom = denominator
        self._prec = 4
        self._expanded = False
        self._reduced = False

        # Only try to eliminate scalar values.
        # TODO: If we divide by a float, we could add the inverse to
        # the numerator as a product, but I don't know if this is
        # efficient since it will involve creating a new object.
        if denominator._prec == 0 and numerator._prec == 0:  # float
            self.num = create_float(numerator.val / denominator.val)
            # Remove denominator, such that it will be excluded when
            # printing.
            self.denom = None

        # Handle zero.
        if self.val == 0.0:
            # Remove denominator, such that it will be excluded when
            # printing
            self.denom = None

        # Compute the representation now, such that we can use it
        # directly in the __eq__ and __ne__ methods (improves
        # performance a bit, but only when objects are cached).
        if self.denom:
            self._repr = "Fraction(%s, %s)" % (self.num._repr, self.denom._repr)
        else:
            self._repr = "Fraction(%s, %s)" % (self.num._repr, create_float(1)._repr)

        # Use repr as hash value.
        self._hash = hash(self._repr)

    # Print functions.
    def __str__(self):
        "Simple string representation which will appear in the generated code."
        if not self.denom:
            return str(self.num)

        # Get string for numerator and denominator.
        num = str(self.num)
        denom = str(self.denom)

        # Group numerator if it is a fraction, otherwise it should be
        # handled already.
        if self.num._prec == 4:  # frac
            num = format["grouping"](num)

        # Group denominator if it is a fraction or product, or if the
        # value is negative.
        # NOTE: This will be removed by the optimisations later before
        # writing any code.
        if self.denom._prec in (2, 4) or self.denom.val < 0.0:  # prod or frac
            denom = format["grouping"](denom)
        return format["div"](num, denom)

    # Binary operators.
    def __add__(self, other):
        "Addition by other objects."
        # Add two fractions if their denominators are equal by creating
        # (expanded) sum of their numerators.
        if other._prec == 4 and self.denom == other.denom:  # frac
            return create_fraction(create_sum([self.num, other.num]).expand(),
                                   self.denom)
        return create_sum([self, other])

    def __sub__(self, other):
        "Subtract other objects."
        # Return a new sum
        if other._prec == 4 and self.denom == other.denom:  # frac
            num = create_sum([self.num, create_product([FloatValue(-1),
                                                        other.num])]).expand()
            return create_fraction(num, self.denom)
        return create_sum([self, create_product([FloatValue(-1), other])])

    def __mul__(self, other):
        "Multiplication by other objects."
        # NOTE: assuming that we get expanded variables.
        # If product will be zero.
        if self.val == 0.0 or other.val == 0.0:
            return create_float(0)
        # Create new expanded numerator and denominator and use '/' to reduce.
        if other._prec != 4:  # frac
            return (self.num * other) / self.denom
        # If we have a fraction, create new numerator and denominator and use
        # '/' to reduce expression.
        return create_product([self.num, other.num]).expand() / create_product([self.denom, other.denom]).expand()

    def __truediv__(self, other):
        "Division by other objects."
        # If division is illegal (this should definitely not happen).
        if other.val == 0.0:
            error("Division by zero.")

        # If fraction will be zero.
        if self.val == 0.0:
            return self.vrs[0]

        # The only thing that we shouldn't need to handle is division by other
        # Fractions
        if other._prec == 4:
            error("Did not expected to divide by fraction.")

        # Handle division by FloatValue, Symbol, Product and Sum in the same
        # way i.e., multiply other by the donominator and use division
        # (__div__ or other) in order to (try to) reduce the expression.
        # TODO: Is it better to first try to divide the numerator by other,
        # if a Fraction is the return value, then multiply the denominator of
        # that value by denominator of self. Otherwise the reduction was
        # successful and we just use the denom of self as denominator.
        return self.num / (other * self.denom)

    __div__ = __truediv__

    # Public functions.
    def expand(self):
        "Expand the fraction expression."

        # If fraction is already expanded, simply return the expansion.
        if self._expanded:
            return self._expanded

        # If we don't have a denominator just return expansion of numerator.
        if not self.denom:
            return self.num.expand()

        # Expand numerator and denominator.
        num = self.num.expand()
        denom = self.denom.expand()

        # TODO: Is it too expensive to call expand in the below?
        # If both the numerator and denominator are fractions, create new
        # numerator and denominator and use division to possibly reduce the
        # expression.
        if num._prec == 4 and denom._prec == 4:  # frac
            new_num = create_product([num.num, denom.denom]).expand()
            new_denom = create_product([num.denom, denom.num]).expand()
            self._expanded = new_num / new_denom
        # If the numerator is a fraction, multiply denominators and use
        # division to reduce expression.
        elif num._prec == 4:  # frac
            new_denom = create_product([num.denom, denom]).expand()
            self._expanded = num.num / new_denom
        # If the denominator is a fraction multiply by the inverse and
        # use division to reduce expression.
        elif denom._prec == 4:  # frac
            new_num = create_product([num, denom.denom]).expand()
            self._expanded = new_num / denom.num
        # Use division to reduce the expression, no need to call expand().
        else:
            self._expanded = num / denom
        return self._expanded

    def get_unique_vars(self, var_type):
        "Get unique variables (Symbols) as a set."
        # Simply get the unique variables from numerator and denominator.
        var = self.num.get_unique_vars(var_type)
        var.update(self.denom.get_unique_vars(var_type))
        return var

    def get_var_occurrences(self):
        """Determine the number of minimum number of times all variables
        occurs in the expression simply by calling the function on the
        numerator.

        """
        return self.num.get_var_occurrences()

    def ops(self):
        "Return number of operations needed to evaluate fraction."
        # If we have a denominator, add the operations and +1 for '/'.
        if self.denom:
            return self.num.ops() + self.denom.ops() + 1
        # Else we just return the number of operations for the
        # numerator.
        return self.num.ops()

    def reduce_ops(self):
        # Try to reduce operations by reducing the numerator and
        # denominator.
        # FIXME: We assume expanded variables here, so any common
        # variables in the numerator and denominator are already
        # removed i.e, there is no risk of encountering (x + x*y) / x
        # -> x*(1 + y)/x -> (1 + y).
        if self._reduced:
            return self._reduced
        num = self.num.reduce_ops()
        # Only return a new Fraction if we still have a denominator.
        if self.denom:
            self._reduced = create_fraction(num, self.denom.reduce_ops())
        else:
            self._reduced = num
        return self._reduced

    def reduce_var(self, var):
        "Reduce the fraction by another variable through division of numerator."
        # We assume that this function is only called by reduce_ops, such that
        # we just need to consider the numerator.
        return create_fraction(self.num / var, self.denom)

    def reduce_vartype(self, var_type):
        """Reduce expression with given var_type. It returns a tuple (found,
        remain), where 'found' is an expression that only has
        variables of type == var_type. If no variables are found,
        found=(). The 'remain' part contains the leftover after
        division by 'found' such that: self = found*remain.

        """

        # Reduce the numerator by the var type.
        if self.num._prec == 3:
            foo = self.num.reduce_vartype(var_type)
            if len(foo) == 1:
                num_found, num_remain = foo[0]
            else:
                # meg: I have only a marginal idea of what I'm doing
                # here!
                new_sum = []
                for num_found, num_remain in foo:
                    if num_found == ():
                        new_sum.append(create_fraction(num_remain, self.denom))
                    else:
                        new_sum.append(create_fraction(create_product([num_found, num_remain]),
                                                       self.denom))
                return create_sum(new_sum).expand().reduce_vartype(var_type)
        else:
            foo = self.num.reduce_vartype(var_type)
            if len(foo) != 1:
                raise RuntimeError("This case is not handled")
            num_found, num_remain = foo[0]

        # If the denominator is not a Sum things are straightforward.
        denom_found = None
        denom_remain = None
        if self.denom._prec != 3:  # sum
            foo = self.denom.reduce_vartype(var_type)
            if len(foo) != 1:
                raise RuntimeError("This case is not handled")
            denom_found, denom_remain = foo[0]

        # If we have a Sum in the denominator, all terms must be
        # reduced by the same terms to make sense
        else:
            remain = []
            for m in self.denom.vrs:
                foo = m.reduce_vartype(var_type)
                d_found, d_remain = foo[0]
                # If we've found a denom, but the new found is
                # different from the one already found, terminate loop
                # since it wouldn't make sense to reduce the fraction.
                # TODO: handle I0/((I0 + I1)/(G0 + G1) + (I1 + I2)/(G1
                # + G2))
                # better than just skipping.
                if len(foo) != 1 or (denom_found is not None and repr(d_found) != repr(denom_found)):
                    # If the denominator of the entire sum has a type
                    # which is lower than or equal to the vartype that
                    # we are currently reducing for, we have to move
                    # it outside the expression as well.
                    # TODO: This is quite application specific, but I
                    # don't see how we can do it differently at the
                    # moment.
                    if self.denom.t <= var_type:
                        if not num_found:
                            num_found = create_float(1)
                        return [(create_fraction(num_found, self.denom),
                                 num_remain)]
                    else:
                        # The remainder is always a fraction
                        return [(num_found, create_fraction(num_remain,
                                                            self.denom))]

                # Update denom found and add remainder.
                denom_found = d_found
                remain.append(d_remain)

            # There is always a non-const remainder if denominator was a sum.
            denom_remain = create_sum(remain)
#        print "den f: ", denom_found
#        print "den r: ", denom_remain
        # If we have found a common denominator, but no found numerator,
        # create a constant.
        # TODO: Add more checks to avoid expansion.
        found = None
        # There is always a remainder.
        remain = create_fraction(num_remain, denom_remain).expand()
#        print "remain: ", repr(remain)

        if num_found:
            if denom_found:
                found = create_fraction(num_found, denom_found)
            else:
                found = num_found
        else:
            if denom_found:
                found = create_fraction(create_float(1), denom_found)
            else:
                found = ()
        return [(found, remain)]
