"Some simple functions for manipulating expressions symbolically"

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2009-07-12 -- 2009-07-15"
__copyright__ = "Copyright (C) 2009 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# FFC common modules
#from ffc.common.log import debug, error

from new_symbol import EPS, CONST, format

def set_format(_format):
    global format
    format = _format

class FloatValue(object):
#    __slots__ = ("val", "t", "_hash")
    def __init__(self, value):
        """Initialise a FloatValue object it contains a:
        val  - float, holds value of object
        t    - Type, always CONST for float."""

        self.val = float(value)
        self.t = CONST

        # TODO: Use cache for hash instead
        self._hash = False

        if abs(value) <  EPS:
            self.val = 0.0

    # Print functions
    def __repr__(self):
        "Representation for debugging"
        return "FloatValue(%s)" % format["floating point"](self.val)

    def __str__(self):
        "Simple string representation"
        return format["floating point"](self.val)

    # Hash (for lookup in {})
    def __hash__(self):
        "Use repr() to compute hash."
        if self._hash:
            return self._hash
        self._hash = hash(repr(self))
        return self._hash

    # Comparison
    def __eq__(self, other):
        "Equal if they are both float values"
        if repr(self) == repr(other):
            return True
        return False

    def __ne__(self, other):
        "Opposite of __eq__"
        return not self == other

    def __lt__(self, other):
        """A float value is always smallest compared to other objects"""
        # If we have two float values, compare numbers
        if isinstance(other, FloatValue):
            return self.val < other.val
        # FloatValue is always lowest
        return True

    def __gt__(self, other):
        "Opposite of __lt__"
        # If we have two float values, compare numbers
        if isinstance(other, FloatValue):
            return self.val > other.val
        # FloatValue is always lowest
        return False

    # Binary operators
    def __add__(self, other):
        "Addition by other objects"
        # NOTE: We expect expanded objects here
        # This is only well-defined if other is a float or if self.val == 0
        if isinstance(other, FloatValue):
            return FloatValue(self.val+other.val)
        elif self.val == 0.0:
            return other
        # Addition is not defined if self is not zero
        # TODO: We could just return a Sum?
        raise RuntimeError("Can only add two floats, or other to zero")

    def __mul__(self, other):
        "Multiplication by other objects"
        # NOTE: We expect expanded objects here i.e., Product([FloatValue])
        # should not be present
        # Only handle case where other is a float, else let the other
        # object handle the multiplication
        if isinstance(other, FloatValue):
            return FloatValue(self.val*other.val)
        return other.__mul__(self)

    def __div__(self, other):
        "Division by other objects"

        # If division is illegal (this should definitely not happen)
        if other.val == 0.0:
            raise RuntimeError("Division by zero")

        # TODO: Should we also support division by fraction for generality?
        # It should not be needed by this module
        if isinstance(other, Fraction):
            raise RuntimeError("Did not expected to divide by fraction")

        # If fraction will be zero
        if self.val == 0.0:
            return self

        # NOTE: We expect expanded objects here i.e., Product([FloatValue])
        # should not be present
        # Handle types appropriately
        if isinstance(other, FloatValue):
            return FloatValue(self.val/other.val)
        # If other is a symbol, return a simple fraction
        elif isinstance(other, Symbol):
            return Fraction(self, other)
        # Don't handle division by sum
        elif isinstance(other, Sum):
            # TODO: Here we could do: 4 / (2*x + 4*y) -> 2/(x + 2*y)
            return Fraction(self, other)

        # If other is a product, remove any float value to avoid
        # 4 / (2*x), this will return 2/x
        val = 1.0
        for v in other.vrs:
            if isinstance(v, FloatValue):
                val *= v.val
        # If we had any floats, create new numerator and only use 'real' variables
        # from the product in the denominator
        if val != 1.0:
            # Check if we need to create a new denominator
            # TODO: Just use other.vrs[1:] instead
            if len(other.get_vrs()) > 1:
                return Fraction(FloatValue(self.val/val), Product(other.get_vrs()))
            # TODO: Because we expect all products to be expanded we shouldn't need
            # to check for this case, just use other.vrs[1]
            elif len(other.get_vrs()) == 1:
                return Fraction(FloatValue(self.val/val), other.vrs[1])
            raise RuntimeError("No variables left in denominator")

        # Nothing left to do
        return Fraction(self, other)

    # Public functions
    def ops(self):
        "Return number of operations to compute the float (always zero)."
        # Just return 0
        return 0

    def expand(self):
        "Expand the expression for a float. (Expanded by construction)."
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
        # A float is not a variable
        return set()

    def reduce_ops(self):
        "Reduce number of operations to evaluate float."
        # Nothing to be done.
        return self

    def get_var_occurrences(self):
        """Determine the number of times all variables occurs in the expression.
        Returns a dictionary of variables and the number of times they occur."""

        # There is only one float value (if it is not -1 or 1)
        if abs(self.val) == 1.0:
            return {}
        return {self:1}

    def reduce_var(self, var):
        "Reduce the float value by another variable through division"
        return self/var

from symbol_obj     import Symbol
from product_obj    import Product
from sum_obj        import Sum
from fraction_obj   import Fraction

