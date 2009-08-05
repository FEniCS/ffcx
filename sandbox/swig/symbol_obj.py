"Some simple functions for manipulating expressions symbolically"

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2009-07-12 -- 2009-07-15"
__copyright__ = "Copyright (C) 2009 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# FFC common modules
#from ffc.common.log import debug, error

from new_symbol import type_to_string

class Symbol(object):
    def __init__(self, variable, symbol_type, base_expr=None, base_op=0):
        """Initialise a Symbols object it contains a:
        val       - float, holds value of object (always 1 for symbol)
        t         - Type, one of CONST, GEO, IP, BASIS
        v         - string, variable name
        base_expr - Other expression type like 'x*y + z'
        base_op   - number of operations for the symbol itself if it's a math
                    operation like std::cos(.) -> base_op = 1."""

        # Dummy value, a symbol is always one
        self.val = 1.0
        self.v = variable
        self.t = symbol_type

        # TODO: Use cache for hash instead
        self._hash = False

        # Needed for symbols like std::cos(x*y + z),
        # where base_expr = x*y + z
        # ops = base_expr.ops() + base_ops = 2 + 1 = 3
        self.base_expr = base_expr
        self.base_op = base_op

        # If type of the base_expr is lower than the given symbol_type change
        # type
        # TODO: Should we raise an error here? Or simply require that one
        # initalise the symbol by Symbol('std::cos(x*y)', (x*y).t, x*y, 1)
        if base_expr and base_expr.t < self.t:
            self.t = base_expr.t

    # Print functions
    def __repr__(self):
        "Representation for debugging"
        if self.base_expr:
            return "Symbol('%s', %s, %s, %d)" % (self.v, type_to_string[self.t], repr(self.base_expr), self.base_op)
        return "Symbol('%s', %s)" % (self.v, type_to_string[self.t])

    def __str__(self):
        "Simple string representation"
        return self.v

    # Hash (for lookup in {})
    def __hash__(self):
        "Use repr as hash"
        if self._hash:
            return self._hash
        # TODO: Will it be OK in all practical cases to use __str__ instead??
        self._hash = hash(repr(self))
        return self._hash

    # Comparison
    def __eq__(self, other):
        "Two symbols are equal if the variable and domain are equal"
        if isinstance(other, Symbol):
            return self.v == other.v and self.t == other.t
        return False

    def __ne__(self, other):
        "Two symbols are not equal if equal is false"
        return not self == other

    def __lt__(self, other):
        "Less than"
        if isinstance(other, FloatValue):
            return False
        # TODO: Just sort by name?
        if isinstance(other, Symbol):
            # First sort by type then by variable name
            if self.t == other.t:
                return self.v < other.v
            elif self.t < other.t:
                return True
            else:
                return False
        # Symbols are always lower than the rest
        return True

    def __gt__(self, other):
        "Greater than"
        if isinstance(other, FloatValue):
            return True
        if isinstance(other, Symbol):
            # First sort by type then by variable name
            if self.t == other.t:
                return self.v > other.v
            elif self.t > other.t:
                return True
            else:
                return False
        # Symbols are always lowest
        return False

    # Binary operators
    def __add__(self, other):
        "Addition of symbols"
        # NOTE: We expect expanded objects and we only expect to add equal
        # symbols, if other is a product, try to let product handle the addition
        # TODO: Should we also support addition by other objects for generality?
        # Add x + x by returning 2*x
        if self == other:
            return Product([FloatValue(2), self])
        elif isinstance(other, Product):
            return other.__add__(self)
        raise RuntimeError("Not implemented")

    def __mul__(self, other):
        "Multiplication of symbols by other objects"

        # NOTE: We assume expanded objects
        # If product will be zero
        if self.val == 0.0 or other.val == 0.0:
            return FloatValue(0)

        # If other is Sum or Fraction let them handle the multiply
        if isinstance(other, (Sum, Fraction)):
            return other.__mul__(self)

        # If other is a float or symbol, create simple product
        if isinstance(other, (FloatValue, Symbol)):
            return Product([self, other])

        # Else add variables from product
        return Product([self] + other.vrs)

    def __div__(self, other):
        "Division of symbols by other objects"

        # NOTE: We assume expanded objects
        # If division is illegal (this should definitely not happen)
        if other.val == 0.0:
            raise RuntimeError("Division by zero")

        # Return 1 if the two symbols are equal
        if self == other:
            return FloatValue(1)

        # If other is a Sum we can only return a fraction
        # TODO: Refine this later such that x / (x + x*y) -> 1 / (1 + y)?
        if isinstance(other, Sum):
            return Fraction(self, other)

        # Handle division by FloatValue, Symbol, Product and Fraction
        # Create numerator and list for denominator
        num = [self]
        denom = []

        # Add floatvalue, symbol and products to the list of denominators
        if isinstance(other, (FloatValue, Symbol)):
            denom.append(other)
        elif isinstance(other, Product):
            denom += other.vrs
        # fraction
        else:
            # TODO: Should we also support division by fraction for generality?
            # It should not be needed by this module
            raise RuntimeError("Did not expected to divide by fraction")

        # Remove one instance of self in numerator and denominator if
        # present in denominator i.e., x/(x*y) --> 1/y
        if self in denom:
            denom.remove(self)
            num.remove(self)

        # Loop entries in denominator and move float value to numerator
        new_denom = []
        append_denom = new_denom.append
        for d in denom:
            # Add the inverse of a float to the numerator and continue
            if isinstance(d, FloatValue):
                num.append(FloatValue(1.0/other.val))
                continue
            append_denom(d)

        # Create appropriate return value depending on remaining data
        # Can only be for x / (2*y*z) -> 0.5*x / (y*z)
        if len(num) > 1:
            num = Product(num)
        # x / (y*z) -> x/(y*z)
        elif num:
            num = num[0]
        # else x / (x*y) -> 1/y
        else:
            num = FloatValue(1)

        # If we have a long denominator, create product and fraction
        if len(new_denom) > 1:
            return Fraction(num, Product(new_denom))
        # If we do have a denominator, but only one variable don't create a
        # product, just return a fraction using the variable as denominator
        elif new_denom:
            return Fraction(num, new_denom[0])
        # If we don't have any donominator left, return the numerator
        # x / 2.0 -> 0.5*x
        return num

    # Public functions
    def ops(self):
        "Returning the number of floating point operation for symbol"

        # Get base ops, typically 1 for sin() and then add the operations
        # for the base (sin(2*x + 1)) --> 2 + 1
        if self.base_expr:
            return self.base_op + self.base_expr.ops()
        return self.base_op

    def expand(self):
        "Expand the expression for a symbol. (Expanded by construction)."
        # Nothing to be done
        return self

    def reduce_vartype(self, var_type):
        """Reduce expression with given var_type. It returns a tuple
        (found, remain), where 'found' is an expression that only has variables
        of type == var_type. If no variables are found, found=(). The 'remain'
        part contains the leftover after division by 'found' such that:
        self = found*remain."""

        if self.t == var_type:
            return (self, FloatValue(1))
        return ((), self)

    def get_unique_vars(self, var_type):
        "Get unique variables (Symbols) as a set"
        # Return self if type matches, also return base expression variables
        s = set()
        if self.t == var_type:
            s.add(self)
        if self.base_expr:
            s.update(self.base_expr.get_unique_vars(var_type))
        return s

    def reduce_ops(self):
        "Reduce number of operations to evaluate symbol."
        # Nothing to be done.
        return self

    def get_var_occurrences(self):
        """Determine the number of times all variables occurs in the expression.
        Returns a dictionary of variables and the number of times they occur."""

        # There is only one symbol
        return {self:1}

    def reduce_var(self, var):
        "Reduce the symbol by another variable through division"
        return self/var

from floatvalue_obj import FloatValue
from product_obj    import Product
from sum_obj        import Sum
from fraction_obj   import Fraction

