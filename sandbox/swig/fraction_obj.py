"This file implements a class to represent a fraction."

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2009-07-12 -- 2009-08-08"
__copyright__ = "Copyright (C) 2009 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# FFC common modules
#from ffc.common.log import error

from new_symbol import create_float, create_product, create_sum, create_fraction
from expr import Expr

#import psyco
#psyco.full()

# TODO: This function is needed to avoid passing around the 'format', but could
# it be done differently?
def set_format(_format):
    global format
    format = _format

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
            raise RuntimeError("Division by zero.")

        # Initialise all variables.
        self.val = numerator.val
        self.t = min([numerator.t, denominator.t])
        self.num = numerator
        self.denom = denominator
        self._prec = 4
        self._expanded = False
        self._reduced = False

        # Only try to eliminate scalar values.
        # TODO: If we divide by a float, we could add the inverse to the
        # numerator as a product, but I don't know if this is efficient
        # since it will involve creating a new object.
        if denominator._prec == 0 and numerator._prec == 0: # float
            self.num = create_float(numerator.val/denominator.val)
            # Remove denominator, such that it will be excluded when printing.
            self.denom = None

        # Handle zero.
        if self.val == 0.0:
            # Remove denominator, such that it will be excluded when printing
            self.denom = None

        # Compute the representation now, such that we can use it directly
        # in the __eq__ and __ne__ methods (improves performance a bit, but
        # only when objects are cached).
        if self.denom:
            self._repr = "Fraction(%s, %s)" %(self.num._repr, self.denom._repr)
        else:
            self._repr = "Fraction(%s, %s)" %(self.num._repr, create_float(1)._repr)

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

        # Group numerator if it is a fraction, otherwise it should be handled already.
        if self.num._prec == 4: # frac
            num = format["grouping"](num)

        # Group denominator if it is a fraction or product, or if the value is negative.
        # NOTE: This will be removed by the optimisations later before writing any code.
        if self.denom._prec in (2, 4) or self.denom.val < 0.0: # prod or frac
            denom = format["grouping"](denom)
        return num + format["division"] + denom

    # Binary operators.
    def __add__(self, other):
        "Addition by other objects."
        # Add two fractions if their denominators are equal by creating
        # (expanded) sum of their numerators.
        if other._prec == 4 and self.denom == other.denom: # frac
            return create_fraction(create_sum([self.num, other.num]).expand(), self.denom)
        else:
            raise RuntimeError("Not implemented.")

    def __mul__(self, other):
        "Multiplication by other objects."
        # NOTE: assuming that we get expanded variables.
        # If product will be zero.
        if self.val == 0.0 or other.val == 0.0:
            return create_float(0)

        # Create new expanded numerator and denominator and use '/' to reduce.
        if other._prec != 4: # frac
            return create_product([self.num, other]).expand()/self.denom
        # If we have a fraction, create new numerator and denominator and use
        # '/' to reduce expression.
        return create_product([self.num, other.num]).expand()/create_product([self.denom, other.denom]).expand()

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
        if num._prec == 4 and denom._prec == 4: # frac
            new_num = create_product([num.num, denom.denom]).expand()
            new_denom = create_product([num.denom, denom.num]).expand()
            self._expanded = new_num/new_denom
        # If the numerator is a fraction, multiply denominators and use
        # division to reduce expression.
        elif num._prec == 4: # frac
            new_denom = create_product([num.denom, denom]).expand()
            self._expanded = num.num/new_denom
        # If the denominator is a fraction multiply by the inverse and
        # use division to reduce expression.
        elif denom._prec == 4: # frac
            new_num = create_product([num, denom.denom]).expand()
            self._expanded = new_num/denom.num
        # Use division to reduce the expression, no need to call expand().
        else:
            self._expanded = num/denom
        return self._expanded

    def get_unique_vars(self, var_type):
        "Get unique variables (Symbols) as a set."
        # Simply get the unique variables from numerator and denominator.
        var = self.num.get_unique_vars(var_type)
        var.update(self.denom.get_unique_vars(var_type))
        return var

    def get_var_occurrences(self):
        """Determine the number of minimum number of times all variables occurs
        in the expression simply by calling the function on the numerator."""
        return self.num.get_var_occurrences()

    def ops(self):
        "Return number of operations needed to evaluate fraction."
        # If we have a denominator, add the operations and +1 for '/'.
        if self.denom:
            return self.num.ops() + self.denom.ops() + 1
        # Else we just return the number of operations for the numerator.
        return self.num.ops()

    def reduce_ops(self):
        # Try to reduce operations by reducing the numerator and denominator.
        # FIXME: We assume expanded variables here, so any common variables in
        # the numerator and denominator are already removed i.e, there is no
        # risk of encountering (x + x*y) / x -> x*(1 + y)/x -> (1 + y).
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
        return create_fraction(self.num/var, self.denom)

    def reduce_vartype(self, var_type):
        """Reduce expression with given var_type. It returns a tuple
        (found, remain), where 'found' is an expression that only has variables
        of type == var_type. If no variables are found, found=(). The 'remain'
        part contains the leftover after division by 'found' such that:
        self = found*remain."""

        # NOTE: We expect self to be expanded at this point.
        # Reduce the numerator by the var type (should be safe, since the
        # expand() should have eliminated all sums in the numerator).
        num_found, num_remain = self.num.reduce_vartype(var_type)

#        # TODO: Remove this test later, expansion should have taken care of
#        # no denominator.
#        if not self.denom:
#            raise RuntimeError("This fraction should have been expanded.")

        # If the denominator is not a Sum things are straightforward.
        denom_found = None
        denom_remain = None
        if self.denom._prec != 3: # sum
            denom_found, denom_remain = self.denom.reduce_vartype(var_type)

        # If we have a Sum in the denominator, all terms must be reduced by
        # the same terms to make sense.
        else:
            remain = []
            for m in self.denom.vrs:
                d_found, d_remain = m.reduce_vartype(var_type)
                # If we've found a denom, but the new found is different from
                # the one already found, terminate loop since it wouldn't make
                # sense to reduce the fraction.
                if denom_found != None and repr(d_found) != repr(denom_found):
                    # In case we did not find any variables of given type in the numerator
                    # declare a constant. We always have a remainder.
                    return (num_found, create_fraction(num_remain, self.denom))

                # Update denom found and add remainder.
                denom_found = d_found
                remain.append(d_remain)

            # There is always a non-const remainder if denominator was a sum.
            denom_remain = create_sum(remain)

        # If we have found a common denominator, but no found numerator,
        # create a constant.
        # TODO: Add more checks to avoid expansion.
        found = None
        # There is always a remainder.
        remain = create_fraction(num_remain, denom_remain).expand()

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
        return (found, remain)

from floatvalue_obj import FloatValue
from symbol_obj     import Symbol
from product_obj    import Product
from sum_obj        import Sum

