"Object to represent a float."

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2009-07-12 -- 2009-08-07"
__copyright__ = "Copyright (C) 2009 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# FFC common modules
#from ffc.common.log import debug, error

from new_symbol import CONST, format, create_float, create_product, create_fraction

import psyco
psyco.full()

def set_format(_format):
    global format
    format = _format
    global format_float
    format_float = format["floating point"]
    global EPS
    EPS = format["epsilon"]

class FloatValue(object):
    __slots__ = ("val", "t", "_class", "_hash", "_repr")
    def __init__(self, value):
        """Initialise a FloatValue object it contains a:
        val     - float, holds value of object.
        t       - Type, always CONST for FloatValue.
        _class  - str, the class type, equal to 'float'.
        _repr   - str, string value of __repr__(), we only compute it once.
        _hash   - int, hash value of __hash__(), we only compute it once."""

        # Initialise value, type and class
        self.val = float(value)
        self.t = CONST
        self._class = "float"

        # Handle zero value explicitly
        if abs(value) <  EPS:
            self.val = 0.0

        # Compute the representation now, such that we can use it directly
        # in the __eq__ and __ne__ methods (improves performance a bit, but
        # only when objects are cached).
        self._repr = "FloatValue(%s)" % format_float(self.val)

        # Use repr as hash value
        self._hash = hash(self._repr)

    # Print functions
    def __repr__(self):
        "Representation for debugging."
        return self._repr

    def __str__(self):
        "Simple string representation."
        return format_float(self.val)

    def __hash__(self):
        "Hash (for lookup in {})."
        return self._hash

    # Comparison
    def __eq__(self, other):
        "==, True if representations are equal."
        # NOTE: 'if other' is not safe since any 'other' which is not one of
        # the symbolic objects will crash this function.
        # However, for speed it is OK to use since anything we get here should
        # be symbolics.
        if other:
            return self._repr == other._repr
        return False

    def __ne__(self, other):
        "!=, True if representations are not equal."
        # NOTE: See above.
        if other:
            return self._repr != other._repr
        return True

    def __lt__(self, other):
        """<, a float value is always smallest compared to other objects.
        If we have two floats compare values."""
        if other._class == "float":
            return self.val < other.val
        return True

    def __gt__(self, other):
        ">, opposite of __lt__."
        if other._class == "float":
            return self.val > other.val
        return False

    # Binary operators
    def __add__(self, other):
        "Addition by other objects."
        # NOTE: We expect expanded objects here
        # This is only well-defined if other is a float or if self.val == 0
        if other._class == "float":
            return create_float(self.val+other.val)
        elif self.val == 0.0:
            return other
        # Addition is not defined if self is not zero
        # TODO: We could just return a Sum?
        raise RuntimeError("Can only add two floats, or other to zero.")

    def __mul__(self, other):
        "Multiplication by other objects."
        # NOTE: We expect expanded objects here i.e., Product([FloatValue])
        # should not be present
        # Only handle case where other is a float, else let the other
        # object handle the multiplication
        if other._class == "float":
            return create_float(self.val*other.val)
        return other.__mul__(self)

    def __div__(self, other):
        "Division by other objects."

        # If division is illegal (this should definitely not happen)
        if other.val == 0.0:
            raise RuntimeError("Division by zero")

        # TODO: Should we also support division by fraction for generality?
        # It should not be needed by this module
        if other._class == "frac":
            raise RuntimeError("Did not expected to divide by fraction")

        # If fraction will be zero
        if self.val == 0.0:
            return self

        # NOTE: We expect expanded objects here i.e., Product([FloatValue])
        # should not be present
        # Handle types appropriately
        if other._class == "float":
            return create_float(self.val/other.val)
        # If other is a symbol, return a simple fraction
        elif other._class == "sym":
            return create_fraction(self, other)
        # Don't handle division by sum
        elif other._class == "sum":
            # TODO: Here we could do: 4 / (2*x + 4*y) -> 2/(x + 2*y)
            return create_fraction(self, other)

        # If other is a product, remove any float value to avoid
        # 4 / (2*x), this will return 2/x
        val = 1.0
        for v in other.vrs:
            if v._class == "float":
                val *= v.val
        # If we had any floats, create new numerator and only use 'real' variables
        # from the product in the denominator
        if val != 1.0:
            # Check if we need to create a new denominator
            # TODO: Just use other.vrs[1:] instead
            if len(other.get_vrs()) > 1:
                return create_fraction(create_float(self.val/val), create_product(other.get_vrs()))
            # TODO: Because we expect all products to be expanded we shouldn't need
            # to check for this case, just use other.vrs[1]
            elif len(other.get_vrs()) == 1:
                return create_fraction(create_float(self.val/val), other.vrs[1])
            raise RuntimeError("No variables left in denominator")

        # Nothing left to do
        return create_fraction(self, other)

    # Public functions
    def ops(self):
        "Return number of operations to compute the float (always zero)."
        # Just return 0
        # NOTE: This means that minus in e.g., -2  and -2*x is not counted.
        return 0

    def expand(self):
        "Expand the expression for a float. (A float is expanded by construction)."
        # Nothing to be done
        return self

    def reduce_vartype(self, var_type):
        """Reduce expression with given var_type. It returns a tuple
        (found, remain), where 'found' is an expression that only has variables
        of type == var_type. If no variables are found, found=(). The 'remain'
        part contains the leftover after division by 'found' such that:
        self = found*remain."""

        if self.t == var_type:
            return (self, create_float(1))
        return ((), self)

    def get_unique_vars(self, var_type):
        "Get unique variables (Symbols) as a set."
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
        "Reduce the float value by another variable through division."
        return self/var

from symbol_obj     import Symbol
from product_obj    import Product
from sum_obj        import Sum
from fraction_obj   import Fraction

