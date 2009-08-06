"Some simple functions for manipulating expressions symbolically"

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2009-07-12 -- 2009-07-15"
__copyright__ = "Copyright (C) 2009 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# FFC common modules
#from ffc.common.log import debug, error

from new_symbol import create_float, create_product, create_sum, create_fraction

#import psyco
#psyco.full()

def set_format(_format):
    global format
    format = _format

class Fraction(object):
    __slots__ = ("val", "t", "num", "denom", "expanded", "reduced", "_class", "_hash", "_repr")
    def __init__(self, numerator, denominator):
        """Initialise a Fraction object, the class contains:
        val   - float, value of fraction (equal to value of numerator)
        num   - expr, the numerator
        denom - expr, the denominator
        t     - Type, one of CONST, GEO, IP, BASIS. It is equal to the lowest
                type of its members"""

        # TODO: Could add 'expanded' to this object too for speed?
        # Check for illegal division
        if denominator.val == 0.0:
            raise RuntimeError("Division by zero")

        # Store the numerator and denominator
        self.num = numerator
        self.denom = denominator

        # Initialise value
        self.val = self.num.val

        # Initialise type
        self.t = min([self.num.t, self.denom.t])

        self.expanded = False
        self.reduced = False

        # Initialise class type
        self._class = "frac"

        # TODO: Use cache for hash instead
        self._hash = False

        self._repr = False
#        self._unique_vars = {}
#        self._reduce_vartype = {}

        # Only try to eliminate scalar values
        # TODO: If we divide by a float, we could add the inverse to the
        # numerator as a product, but I don't know if this is efficient
        # since it will involve creating a new object
#        if isinstance(denominator, FloatValue) and isinstance(numerator, FloatValue):
        if denominator._class == "float" and numerator._class == "float":
#            self.num = FloatValue(numerator.val/denominator.val)
            self.num = create_float(numerator.val/denominator.val)
            # Remove denominator, such that it will be excluded when printing
            self.denom = None

        # Handle zero
        if self.val == 0.0:
            # Remove denominator, such that it will be excluded when printing
            self.denom = None

    # Print functions
    def __repr__(self):
        "Representation for debugging"
        if not self._repr:
            if self.denom:
                self._repr = "Fraction(%s, %s)" %(repr(self.num), repr(self.denom))
            else:
#                self._repr = "Fraction(%s, %s)" %(repr(self.num), repr(FloatValue(1)))
                self._repr = "Fraction(%s, %s)" %(repr(self.num), repr(create_float(1)))
        return self._repr


    def __str__(self):
        "Simple string representation"
        if not self.denom:
            return str(self.num)

        # Get string for numerator and denominator
        num = str(self.num)
        denom = str(self.denom)

        # Group numerator if it is a fraction, otherwise it should be handled
        # already
#        if isinstance(self.num, Fraction):
        if self.num._class == "frac":
            num = format["grouping"](num)

        # Group denominator if it is a fraction or product, or if the value is
        # negative.
        # NOTE: This will be removed by the optimisations later before writing
        # any code
#        if isinstance(self.denom, (Product, Fraction)) or self.denom.val < 0.0:
        if self.denom._class in ("prod", "frac") or self.denom.val < 0.0:
            denom = format["grouping"](denom)

        return num + format["division"] + denom

    # Hash (for lookup in {})
    def __hash__(self):
        "Use repr as hash"
        if self._hash:
            return self._hash
        self._hash = hash(repr(self))
        return self._hash

    def __eq__(self, other):
        # Fractions are equal if their denominator and numerator are equal
#        return repr(self) == repr(other)
#        if isinstance(other, Fraction):
        if other and other._class == "frac":
            return self.denom == other.denom and self.num == other.num
        return False

    def __ne__(self, other):
        "Two fractions are not equal if equal is false"
#        return repr(self) != repr(other)
        return not self == other

    def __lt__(self, other):
        # Fractions are always greater than
#        if isinstance(other, Fraction):
        if other._class == "frac":
            if self.num < other.num:
                return True
            elif self.num == other.num and self.denom < other.denom:
                return True
        return False

    def __gt__(self, other):
#        if isinstance(other, Fraction):
        if other._class == "frac":
            if self.num > other.num:
                return True
            elif self.num == other.num and self.denom > other.denom:
                return True
            else:
                return False
        return True

    # Binary operators
    def __add__(self, other):

        # Add two fractions if their denominators are equal by creating
        # (expanded) sum of their numerators
#        if isinstance(other, Fraction) and self.denom == other.denom:
        if other._class == "frac" and self.denom == other.denom:
#            return Fraction(Sum([self.num, other.num]).expand(), self.denom)
#            return Fraction(create_sum([self.num, other.num]).expand(), self.denom)
            return create_fraction(create_sum([self.num, other.num]).expand(), self.denom)
        else:
            raise RuntimeError("Not implemented")

    def __mul__(self, other):

        # NOTE: assuming that we get expanded variables
        # If product will be zero
        if self.val == 0.0 or other.val == 0.0:
#            return FloatValue(0)
            return create_float(0)

        # Create new expanded numerator and denominator and use '/' to reduce.
#        if not isinstance(other, Fraction):
        if other._class != "frac":
#            return Product([self.num, other]).expand()/self.denom
            return create_product([self.num, other]).expand()/self.denom
        # If we have a fraction, create new numerator and denominator and use
        # '/' to reduce expression
#        return Product([self.num, other.num]).expand()/Product([self.denom, other.denom]).expand()
        return create_product([self.num, other.num]).expand()/create_product([self.denom, other.denom]).expand()

    # Public functions
    def ops(self):
        "Return number of operations needed to evaluate fraction"
        # If we have a denominator, add the operations and +1 for '/'
        if self.denom:
            return self.num.ops() + self.denom.ops() + 1

        # Else we just return the number of operations for the numerator
        return self.num.ops()

    def expand(self):
        "Expand the fraction expression"

        if self.expanded:
            return self.expanded

        # If we don't have a denominator just return expansion of numerator
        if not self.denom:
            return self.num.expand()

        # Expand numerator and denominator
        num = self.num.expand()
        denom = self.denom.expand()

        # TODO: Is it too expensive to call expand in the below?
        # If both the numerator and denominator are fractions, create new
        # numerator and denominator and use division to possibly reduce the
        # expression
#        if isinstance(num, Fraction) and isinstance(denom, Fraction):
        if num._class == "frac" and denom._class == "frac":
#            new_num = Product([num.num, denom.denom]).expand()
            new_num = create_product([num.num, denom.denom]).expand()
#            new_denom = Product([num.denom, denom.num]).expand()
            new_denom = create_product([num.denom, denom.num]).expand()
            self.expanded = new_num/new_denom

        # If the numerator is a fraction, multiply denominators and use
        # division to reduce expression
#        elif isinstance(num, Fraction):
        elif num._class == "frac":
#            new_denom = Product([num.denom, denom]).expand()
            new_denom = create_product([num.denom, denom]).expand()
            self.expanded = num.num/new_denom

        # If the denominator is a fraction multiply by the inverse and
        # use division to reduce expression
#        elif isinstance(denom, Fraction):
        elif denom._class == "frac":
#            new_num = Product([num, denom.denom]).expand()
            new_num = create_product([num, denom.denom]).expand()
            self.expanded = new_num/denom.num

        # Use division to reduce the expression, no need to call expand()
        else:
            self.expanded = num/denom
        return self.expanded


    def reduce_vartype(self, var_type):
        """Reduce expression with given var_type. It returns a tuple
        (found, remain), where 'found' is an expression that only has variables
        of type == var_type. If no variables are found, found=(). The 'remain'
        part contains the leftover after division by 'found' such that:
        self = found*remain."""

#        if var_type in self._reduce_vartype:
#            return self._reduce_vartype[var_type]

        # NOTE: We expect self to be expanded at this point
        # Reduce the numerator by the var type (should be safe, since the
        # expand() should have eliminated all sums in the numerator)
        num_found, num_remain = self.num.reduce_vartype(var_type)

        # TODO: Remove this test later, expansion should have taken care of
        # no denominator
        if not self.denom:
            raise RuntimeError("This fraction should have been expanded")

        # If the denominator is not a Sum things are straightforward
        denom_found = None
        denom_remain = None
#        if not isinstance(self.denom, Sum):
        if self.denom._class != "sum":
            denom_found, denom_remain = self.denom.reduce_vartype(var_type)

        # If we have a Sum in the denominator, all terms must be reduced by
        # the same terms to make sense
        else:
            remain = []
            for m in self.denom.pos + self.denom.neg:
                d_found, d_remain = m.reduce_vartype(var_type)
                # If we've found a denom, but the new found is different from
                # the one already found, terminate loop since it wouldn't make
                # sense to reduce the fraction
                if denom_found != None and str(d_found) != str(denom_found):
                    # In case we did not find any variables of given type in the numerator
                    # declare a constant. We always have a remainder.
#                    return (num_found, Fraction(num_remain, self.denom))
                    return (num_found, create_fraction(num_remain, self.denom))
#                    self._reduce_vartype[var_type] = (num_found, create_fraction(num_remain, self.denom))
#                    return self._reduce_vartype[var_type]

                denom_found = d_found
                remain.append(d_remain)

            # There is always a non-const remainder if denominator was a sum
#            denom_remain = Sum(remain)
            denom_remain = create_sum(remain)

        # If we have found a common denominator, but no found numerator,
        # create a constant
        # TODO: Add more checks to avoid expansion
        found = None
        # There is always a remainder
#        remain = Fraction(num_remain, denom_remain).expand()
        remain = create_fraction(num_remain, denom_remain).expand()

        if num_found:
            if denom_found:
#                found = Fraction(num_found, denom_found)
                found = create_fraction(num_found, denom_found)
            else:
                found = num_found
        else:
            if denom_found:
#                found = Fraction(FloatValue(1), denom_found)
#                found = Fraction(create_float(1), denom_found)
                found = create_fraction(create_float(1), denom_found)
            else:
                found = ()
        return (found, remain)
#        self._reduce_vartype[var_type] = (found, remain)
#        return (found, remain)

    def get_unique_vars(self, var_type):
        "Get unique variables (Symbols) as a set"
#        if var_type in self._unique_vars:
#            return self._unique_vars[var_type]
        var = set()
        # Simply get the unique variables from numerator and denominator
        var.update(self.num.get_unique_vars(var_type))
        var.update(self.denom.get_unique_vars(var_type))
#        self._unique_vars[var_type] = var
        return var

    def reduce_ops(self):
        # Try to reduce operations by reducing the numerator and denominator.
        # FIXME: We assume expanded variables here, so any common variables in
        # the numerator and denominator are already removed i.e, there is no
        # risk of encountering (x + x*y) / x -> x*(1 + y)/x -> (1 + y)
        if self.reduced:
            return self.reduced
        num = self.num.reduce_ops()
        # Only return a new Fraction if we still have a denominator
        if self.denom:
            denom = self.denom.reduce_ops()
#            return Fraction(num, denom)
            self.reduced = create_fraction(num, denom)
        else:
            self.reduced = num
        return self.reduced

    def get_var_occurrences(self):
        """Determine the number of minimum number of times all variables occurs
        in the expression simply by calling the function on the numerator"""
        return self.num.get_var_occurrences()

    def reduce_var(self, var):
        "Reduce the fraction by another variable through division of numerator."
        # We assume that this function is only called by reduce_ops, such that
        # we just need to consider the numerator.
#        return Fraction(self.num/var, self.denom)
        return create_fraction(self.num/var, self.denom)

from floatvalue_obj import FloatValue
from symbol_obj     import Symbol
from product_obj    import Product
from sum_obj        import Sum

