"Some simple functions for manipulating expressions symbolically"

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2009-07-12 -- 2009-07-15"
__copyright__ = "Copyright (C) 2009 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# FFC common modules
#from ffc.common.log import debug, error
#from copy import deepcopy

def set_format(_format):
    global format
    format = _format

class Product(object):
    def __init__(self, variables):
        """Initialise a Product object, the class contains:
        val       - float, holds the value of the object
        t         - Type, one of CONST, GEO, IP, BASIS. It is equal to the lowest
                    type of its members
        vrs       - a list of variables
        expanded  - bool, flag depending whether or not the product needs expansion"""

        self.val = 1.0
        self.vrs = []
        self.expanded = True

        # TODO: Use cache for hash instead
        self._hash = False

        if variables:
            # Remove nested Products and test for expansion
            new_vars = []
            for v in variables:
                # If any value is zero the entire term is zero
                if v.val == 0.0:
                    self.val = 0.0
                    self.vrs = [FloatValue(0.0)]
                    new_vars = []
                    break

                # If we have sums or fractions in the variables the product is
                # not expanded
                if isinstance(v, (Sum, Fraction)):
                    self.expanded = False

                # Take care of product such that we don't create nested products
                if isinstance(v, Product):
                    # If other product is not expanded, we must expand this product later
                    if not v.expanded:
                        self.expanded = False
                    new_vars += v.vrs
                    continue
                new_vars.append(v)

            # Loop variables and collect all floats in one variable
            float_val = 1.0
            for v in new_vars:
                if isinstance(v, FloatValue):
                    float_val *= v.val
                    continue
                self.vrs.append(v)

            # If value is 1 there is no need to include it, unless it is the
            # only parameter left i.e., 2*0.5 = 1
            if float_val and float_val != 1.0:
                fv = FloatValue(float_val)
                self.val = float_val
                self.vrs.append(fv)
            elif float_val == 1.0 and not self.vrs:
                fv = FloatValue(float_val)
                self.val = float_val
                self.vrs.append(fv)

        # If we don't have any variables the product is zero
        else:
            self.vrs = [FloatValue(0)]

        # The type is equal to the lowest variable type
        self.t = min([v.t for v in self.vrs])

        # Sort the variables such that comparisons work
        self.vrs.sort()

    # Print functions
    def __repr__(self):
        "Representation for debugging"
        return "Product([%s])" % ", ".join([repr(v) for v in self.vrs])

    def __str__(self):
        "Simple string representation"

        # If the first float is -1 exlude the 1.
        if isinstance(self.vrs[0], FloatValue) and self.vrs[0].val == -1.0:
            # Join string representation of members by multiplication
            neg = format["subtract"](["",""]).split()[0]
            return neg + format["multiply"]([str(v) for v in self.vrs[1:]])
        return format["multiply"]([str(v) for v in self.vrs])

    # Hash (for lookup in {} and [])
    def __hash__(self):
        "Use repr as hash"
        if self._hash:
            return self._hash
        # TODO: Better way than use hash on all members?
        # NOTE: Exclude FloatValue when computing hash
#        self._hash = hash(tuple([hash(v) for v in self.get_vrs()]))
        self._hash = hash(repr(self))
        return self._hash

    # Comparison
    # TODO: Beware that due to definition of
    # FloatValue and Symbol definitions 3*x == 2*x evaluates to True
    def __eq__(self, other):
        "Two products are equal if their list of variables are equal"
        if isinstance(other, Product):
#            return self.get_vrs() == other.get_vrs()
            return self.vrs == other.vrs
        else:
            return False

    def __ne__(self, other):
        "Two products are not equal if equal is false"
        return not self == other

    def __lt__(self, other):
        # FloatValue and Symbols are always less
        if isinstance(other, (FloatValue, Symbol)):
            return False
        # Compare list of symbols for two products
        elif isinstance(other, Product):
            return self.vrs < other.vrs
        # Products are less than sum and fraction
        return True

    def __gt__(self, other):
        # Symbols are always less
        if isinstance(other, (FloatValue, Symbol)):
            return True
        # Compare list of symbols for two products
        elif isinstance(other, Product):
            return self.vrs > other.vrs
        # Products are less than sum and fraction
        return False

    # Binary operators
    def __add__(self, other):
        # NOTE: Assuming expanded variables
        # If two products are equal, add their float values
        if isinstance(other, Product) and self.get_vrs() == other.get_vrs():
            # Return expanded product, to get rid of 3*x + -2*x -> x, not 1*x
            return Product([FloatValue(self.val + other.val)] + list(self.get_vrs())).expand()
        # if self == 2*x and other == x return 3*x
        elif isinstance(other, Symbol):
            if self.get_vrs() == (other,):
                # Return expanded product, to get rid of -x + x -> 0, not product(0)
                return Product([FloatValue(self.val + 1.0), other]).expand()
            else:
                raise RuntimeError("Not implemented")
        else:
            raise RuntimeError("Not implemented")

    def __mul__(self, other):

        # If product will be zero
        if self.val == 0.0 or other.val == 0.0:
            return FloatValue(0)

        # If other is a Sum or Fraction let them handle it
        if isinstance(other, (Sum, Fraction)):
            return other.__mul__(self)

        # Get copy of variables
        new_prod = [v for v in self.vrs]

        # NOTE: We expect expanded sub-expressions with no nested operators
        # If other is a float or symbol, add to list
        if isinstance(other, (FloatValue, Symbol)):
            new_prod.append(other)
        # Add variables of other product
        else:
            new_prod += other.vrs
        return Product(new_prod)

    def __div__(self, other):

        # If division is illegal (this should definitely not happen)
        if other.val == 0.0:
            raise RuntimeError("Division by zero")

        # If fraction will be zero
        if self.val == 0.0:
            return self.vrs[0]

        # If other is a Sum we can only return a fraction
        # NOTE: Expect that other is expanded i.e., x + x -> 2*x which can be handled
        # TODO: Fix x / (x + x*y) -> 1 / (1 + y)
        # Or should this be handled when reducing a fraction?
        if isinstance(other, Sum):
            return Fraction(self, other)

        # Handle division by FloatValue, Symbol, Product and Fraction
        # NOTE: assuming that we get expanded variables

        # Copy numerator, and create list for denominator
        num = [v for v in self.vrs]
        denom = []

        # Add floatvalue, symbol and products to the list of denominators
        if isinstance(other, (FloatValue, Symbol)):
            denom.append(other)
        elif isinstance(other, Product):
            denom += other.vrs
        # fraction
        else:
            raise RuntimeError("Did not expected to divide by fraction")

        # Loop entries in denominator and remove from numerator (and denominator)
        new_denom = []
        for d in denom:
            # Add the inverse of a float to the numerator and continue
            if isinstance(d, FloatValue):
                num.append(FloatValue(1.0/d.val))
                continue
            if d in num:
                num.remove(d)
            else:
                new_denom.append(d)

        # Create appropriate return value depending on remaining data
        if len(num) > 1:
            # TODO: Make this more efficient?
            # Create product and expand to reduce
            # Product([5, 0.2]) == Product([1]) -> Float(1)
            num = Product(num).expand()
        elif num:
            num = num[0]
        # If all variables in the numerator has been eliminated we only have
        # 1 left (We should then have no denominator left)
        else:
            if new_denom:
                raise RuntimeError("Where did the numbers go?")
            num = FloatValue(1)

        if len(new_denom) > 1:
            new_denom = Product(new_denom)
            return Fraction(num, new_denom)
        elif new_denom:
            return Fraction(num, new_denom[0])
        return num

    # Public functions
    def ops(self):
        "Get the number of operations to compute product"

        # It takes n-1 operations ('*') for a product of n members
        op = len(self.vrs) - 1

        # Loop members and add their count
        for v in self.vrs:
            op += v.ops()

        # Subtract 1, if the first member is -1 i.e., -1*x*y -> x*y is only 1 op.
        if isinstance(self.vrs[0], FloatValue) and self.vrs[0].val == -1.0:
            return op - 1
        return op

    def expand(self):
        "Expand all members of the product"

        # If we just have one variable, return the expansion of it
        # (it is not a Product, so it should be safe). We need this to get
        # rid of Product([Symbol]) type expressions
        if len(self.vrs) == 1:
            return self.vrs[0].expand()

        # If product is already expanded, return self
        if self.expanded:
            return self

        # Sort variables in FloatValue and Symbols and the rest such that
        # we don't call the '*' operator more than we have to
        float_syms = []
        rest = []
        append_syms = float_syms.append
        append_rest = rest.append
        for v in self.vrs:
            if isinstance(v, (FloatValue, Symbol)):
                append_syms(v)
                continue

            append_rest(v.expand())

        # If we have floats or symbols add the symbols to the rest using
        # appropriate object
        # single product (for speed)
        if len(float_syms) > 1:
            append_rest( Product(float_syms) )
        elif float_syms:
            append_rest(float_syms[0])

        # Use __mult__ to reduce list to one single variable
        # TODO: Can this be done more efficiently without creating all the
        # intermediate variables?
        return reduce(lambda x,y: x*y, rest)

    def get_vrs(self):
        "Return all 'real' variables"
        # A product should only have one float value
        # TODO: Use this knowledge directly in other classes
        if isinstance(self.vrs[0], FloatValue):
            return tuple(self.vrs[1:])
        return tuple(self.vrs)

    def reduce_vartype(self, var_type):
        """Reduce expression with given var_type. It returns a tuple
        (found, remain), where 'found' is an expression that only has variables
        of type == var_type. If no variables are found, found=(). The 'remain'
        part contains the leftover after division by 'found' such that:
        self = found*remain."""

        # Sort variables according to type
        found = []
        remains = []
        for v in self.vrs:
            if v.t == var_type:
                found.append(v)
                continue
            remains.append(v)

        # Create appropriate object for found
        if len(found) > 1:
            found = Product(found)
        elif found:
            found = found.pop()
        # We did not find any variables
        else:
            return ((), self)

        # Create appropriate object for remains
        if len(remains) > 1:
            remains = Product(remains)
        elif remains:
            remains = remains.pop()
        # We don't have anything left
        else:
            return (self, FloatValue(1))

        # Return whatever we found
        return (found, remains)

    def get_unique_vars(self, var_type):
        "Get unique variables (Symbols) as a set"
        # Loop all members and get their types
        var = set()
        update = var.update
        for v in self.vrs:
            update(v.get_unique_vars(var_type))
        return var

    def reduce_ops(self):
        "Reduce the number of operations to evaluate the product"
        # It's not possible to reduce a product if it is already expanded and
        # it should be at this stage
        # TODO: Is it safe to return self.expand().reduce_ops() if product is
        # not expanded? And do we want to?
        if self.expanded:
            return self
        raise RuntimeError("Product must be expanded first before we can reduce the number of operations")

    def get_var_occurrences(self):
        """Determine the number of times all variables occurs in the expression.
        Returns a dictionary of variables and the number of times they occur."""

        # TODO: The product should be expanded at this stage, should we check
        # this?

        # Create dictionary and count number of occurrences of each variable
        d = {}
        for v in self.vrs:
            if v in d:
                d[v] += 1
            else:
                d[v] = 1
        return d


#    def num_var(self, var):
#        # The number of variables with a given name is just the sum of all
#        # occurrences of that symbol in the product
#        return sum([v == var for v in self.vrs])


#    def reduce_var(self, var):
#        # Reduce by another variable by division
#        return self/var

#    def get_all_vars(self):
#        return self.vrs


from floatvalue_obj import FloatValue
from symbol_obj     import Symbol
from sum_obj        import Sum
from fraction_obj   import Fraction


