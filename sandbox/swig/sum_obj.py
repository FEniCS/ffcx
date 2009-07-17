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

class Sum(object):
    def __init__(self, variables):
        """Initialise a Sum object, the class contains:
        val - float, numeric value of the object (defualt is 1, but can be
              zero if sum does not contain any members)
        t   - Type, one of CONST, GEO, IP, BASIS. It is equal to the lowest
              type of its members
        pos - list, all positive variables
        neg - list, all negative variables
        """

        # TODO: Could add 'expanded' to this object too for speed?
        self.val = 1.0
        self.pos = []
        self.neg = []

        # TODO: Use cache for hash instead
        self._hash = False

        if variables:
            # Remove nested Sums
            new_vars = []
            for v in variables:
                if isinstance(v, Sum):
                    new_vars += v.pos + v.neg
                    continue
                new_vars.append(v)

            floats = []
            # Loop variables and sort in positive and negative
            # Also collect all floats in 1 variable
            # We don't collect [x, x, x] into 3*x to avoid creating objects
            # Instead we do this when expanding the object
            for v in new_vars:
                # Skip zero terms
                if v.val ==  0.0:
                    continue
                if isinstance(v, FloatValue):
                    floats.append(v)
                elif v.val < 0.0:
                    self.neg.append(v)
                else:
                    self.pos.append(v)

            # Only create new float if we have more than one. Also ignore it
            # if 1 + 1 - 2 = 0
            if len(floats) > 1:
                val = sum([v.val for v in floats])
                if val and val < 0.0:
                    self.neg.append(FloatValue(val))
                elif val:
                    self.pos.append(FloatValue(val))
            elif floats:
                if floats[0].val and floats[0].val < 0.0:
                    self.neg.append(floats[0])
                elif floats[0].val:
                    self.pos.append(floats[0])

        # If we don't have any variables the sum is zero
        else:
            self.pos = [FloatValue(0)]

        # Handle zero value
        if not self.pos + self.neg:
            self.val = 0.0
            self.pos = [FloatValue(0)]

        # Type is equal to the smallest type in both lists
        self.t = min([v.t for v in self.neg + self.pos])

        # Sort variables in positive and negative, (for representation)
        self.pos.sort()
        self.neg.sort()

    # Print functions
    def __repr__(self):
        "Representation for debugging"
        return "Sum([%s])" % ", ".join([repr(v) for v in self.pos + self.neg])

    def __str__(self):
        "Simple string representation"

        # First add all the positive variables using plus, then add all
        # negative variables
        s = format["add"]([str(v) for v in self.pos]) + \
            "".join([str(v) for v in self.neg])

        # Group only if we have more that one variable
        if len(self.pos + self.neg) > 1:
            return format["grouping"](s)
        return s

    # Hash (for lookup in {} and [])
    def __hash__(self):
        "Use repr as hash"
        if self._hash:
            return self._hash
        # TODO: Better way than use hash on all members?
        # NOTE: Exclude FloatValue when computing hash
#        self._hash = hash(tuple([hash(v) for v in self.pos + self.neg]))
        self._hash = hash(repr(self))
        return self._hash

    # Comparison
    # TODO: Beware that due to definition of
    # FloatValue and Symbol definitions 3*x == 2*x evaluates to True
    def __eq__(self, other):
        "Two sums are equal if their list of variables are equal including sign"
        if isinstance(other, Sum):
            return self.pos == other.pos and self.neg == other.neg
        else:
            return False

    def __ne__(self, other):
        "Two sums are not equal if equal is false"
        return not self == other

    def __lt__(self, other):
        # Symbols and products are always less
        if isinstance(other, (FloatValue, Symbol, Product)):
            return False
        # Compare representation, to get members and sign correct
        elif isinstance(other, Sum):
            return self.pos + self.neg < other.pos + other.neg
        # Fractions are always greater
        return True

    def __gt__(self, other):
        # Symbols and products are always less
        if isinstance(other, (FloatValue, Symbol, Product)):
            return True
        # Compare representation, to get members and sign correct
        elif isinstance(other, Sum):
            return self.pos + self.neg > other.pos + other.neg
        return False

    # Binary operators
    def __mul__(self, other):

        # If product will be zero
        if self.val == 0.0 or other.val == 0.0:
            return FloatValue(0)

        # NOTE: We expect expanded sub-expressions with no nested operators
        # Create list of new products using the '*' operator
        # TODO: Is this efficient?
        new_prods = [v*other for v in self.pos + self.neg]

        # Remove zero valued terms
        # TODO: Can this still happen?
        new_prods = [v for v in new_prods if v.val != 0.0]

        # Create new sum
        if not new_prods:
            return FloatValue(0)
        elif len(new_prods) > 1:
            # Expand sum to collect terms
            return Sum(new_prods).expand()
        # TODO: Is it necessary to call expand?      
        return new_prods[0].expand()

    def __div__(self, other):

        # If division is illegal (this should definitely not happen)
        if other.val == 0.0:
            raise RuntimeError("Division by zero")

        # If fraction will be zero
        if self.val == 0.0:
            return FloatValue(0)

        # NOTE: assuming that we get expanded variables
        # If other is a Sum we can only return a fraction
        # TODO: We could check for equal sums if Sum.__eq__ could be trusted
        # As it is now (2*x + y) == (3*x + y), which works for the other things I do
        # NOTE: Expect that other is expanded i.e., x + x -> 2*x which can be handled
        # TODO: Fix (1 + y) / (x + x*y) -> 1 / x
        # Will this be handled when reducing operations on a fraction?
        if isinstance(other, Sum):
            return Fraction(self, other)


        # NOTE: We expect expanded sub-expressions with no nested operators
        # Create list of new products using the '*' operator
        # TODO: Is this efficient?
        new_fracs = [v/other for v in self.pos + self.neg]

        # Remove zero valued terms
        # TODO: Can this still happen?
        new_fracs = [v for v in new_fracs if v.val != 0.0]

        # Create new sum
        # TODO: No need to call expand here, using the '/' operator should have
        # taken care of this
        if not new_fracs:
            return FloatValue(0)
        elif len(new_fracs) > 1:
            return Sum(new_fracs)
        return new_fracs[0]

    # Public functions
    def ops(self):
        "Return number of operations to compute value of sum"
        op = 0
        # Add the number of operations from sub-expressions
        for v in self.pos + self.neg:
            #  +1 for the +/- symbol
            op += v.ops() + 1

        # Subtract one operation as it only takes n-1 ops to sum n members
        return op - 1

    def expand(self):
        "Expand all members of the sum"

        # TODO: This function might need some optimisation
        # Get all variables and expand (is OK because we don't have any nested Sums)
        variables = [v.expand() for v in self.pos + self.neg]

        # Remove nested sums that got created by expanding variables
        new_vars = []
        new_append = new_vars.append
        for v in variables:
            if isinstance(v, Sum):
                new_vars += v.pos + v.neg
                continue
            new_append(v)

        # Sort variables into symbols, products and fractions (add floats
        # directly to new list, will be handled later). Add fractions if
        # possible else add to list.
        new_variables = []
        append = new_variables.append
        syms = []
        sym_append = syms.append
        prods = []
        prod_append = prods.append
        frac_groups = {}
        # TODO: Rather than using '+', would it be more efficient to collect
        # the terms first?
        for v in new_vars:
#            print "v: ", repr(v)
#            if isinstance(v, FloatValue):
            # TODO: Should we also group fractions, or put this in a separate function?
            if isinstance(v, (FloatValue, Fraction)):
                append(v)
            elif isinstance(v, Symbol):
                sym_append(v)
            elif isinstance(v, Product):
                prod_append(v)
            # TODO: put this in another function, cannot group fractions
            # before we have reduced fractions with respect to Types etc.
#            else:
#                # Try to group fractions with similar denominators
#                if v.denom in frac_groups:
#                    frac_groups[v.denom] += v
#                else:
#                    frac_groups[v.denom] = v

        # Sort all variables in groups: [2*x, -7*x], [(x + y), (2*x + 4*y)] etc.
        # First handle product in order to add symbols if possible
        prod_groups = {}
        for v in prods:
#            print "v: ", v
#            print "getvrs: ", v.get_vrs()
            if v.get_vrs() in prod_groups:
                prod_groups[v.get_vrs()] += v
            else:
                prod_groups[v.get_vrs()] = v

        sym_groups = {}
        # Loop symbols and add to appropriate groups
        for v in syms:
            # First try to add to a product group
            if (v,) in prod_groups:
                prod_groups[(v,)] = v
            # Then to other symbols
            elif v in sym_groups:
                sym_groups[v] += v
            # Create a new entry in the symbols group
            else:
                sym_groups[v] = v

        # Loop groups and add to new variable list
        for k,v in sym_groups.items():
            append(v)
        for k,v in prod_groups.items():
            append(v)
#        for k,v in frac_groups.items():
#            append(v)

        if len(new_variables) > 1:
            # Return new sum (will remove multiple instances of floats during construction)
            return Sum(new_variables)
        elif new_variables:
            # If we just have one variable left, return it since it is already expanded
            return new_variables[0]
        raise RuntimeError("Where did the variables go?")

    def reduce_vartype(self, var_type):
        """Reduce expression with given var_type. It returns a list of tuples
        [(found, remain)], where 'found' is an expression that only has variables
        of type == var_type. If no variables are found, found=(). The 'remain'
        part contains the leftover after division by 'found' such that:
        self = Sum([f*r for f,r in self.reduce_vartype(Type)])."""

        found = {}
        # Loop members and reduce them by vartype
        for v in self.pos + self.neg:
            f, r = v.reduce_vartype(var_type)
            if f in found:
                found[f].append(r)
            else:
                found[f] = [r]

        # Create the return value
        returns = []
        append = returns.append
        for f, r in found.items():
            if len(r) > 1:
                # Use expand to group expressions
                r = Sum(r).expand()
            elif r:
                r = r.pop()
            append((f, r))
        return returns

    def get_unique_vars(self, var_type):
        "Get unique variables (Symbols) as a set"
        # Loop all variables of self update the set
        var = set()
        update = var.update
        for v in self.pos + self.neg:
            update(v.get_unique_vars(var_type))
        return var

    def reduce_ops(self):
        "Reduce the number of operations needed to evaluate the sum"

        # NOTE: Assuming that sum has already been expanded
        # TODO: Add test for this and handle case if it is not

        # TODO: This looks expensive, can it be optimised?

        # Loop all variables of the sum and collect the number of common
        # variables that can be factored out.
        common_vars = {}
        for var in self.pos + self.neg:

            # Get dictonary of occurrences
            # TODO: Don't save d, group with following line
            d = var.get_var_occurrences()

            # Add the variable and the number of occurrences to common dictionary
            for k, v in d.items():
                if k in common_vars:
                    common_vars[k].append((v, var))
                else:
                    common_vars[k] = [(v, var)]

        print "\n\nvars: ", common_vars

        # Determine the maximum reduction for each variable
        reduction_vars = []
        do_reduce_terms = []
        max_red = 0
        for k,v in common_vars.items():
            print k, v
            # If the number of expressions that can be reduced is only one
            # there is nothing to be done
            if len(v) > 1:
                # Compute possible reduction
                reduc = len(v)
                occur = [t[0] for t in v]
                min_occur = min(occur)
                max_occur = max(occur)

                if reduc > max_red:
                    max_red = reduc
                    reduction_vars = [k]
                    do_reduce_terms = v
                # If this reducing by this variable reduces the number of
                # operations similarly, check if we can reduce all at once.
                elif reduc == max_red:
                    # TODO: Do we need to do some kind or sorting here?
                    if v == do_reduce_terms:
                        reduction_vars.append(k)
            else:
                continue

        if reduction_vars:
            print "reduction vars: ", reduction_vars

            # Extract terms from tuples
            do_reduce_terms = [v[1] for v in do_reduce_terms]
            print "do_reduce_terms: ", do_reduce_terms

            # Create variable that we will use to reduce the terms
            reduction_var = None
            if len(reduction_vars) > 1:
                reduction_var = Product(reduction_vars)
            else:
                reduction_var = reduction_vars[0]

            # Reduce all terms that need to be reduced
            reduced_terms = [v/reduction_var for v in do_reduce_terms]

            # Create reduced expression
            reduced_expr = None
            if len(reduced_terms) > 1:
                # Try to reduce the reduced terms further
                reduced_expr = Product([reduction_var, Sum(reduced_terms).reduce_ops()])
            else:
                reduced_expr = Product(reduction_var, reduced_terms[0])

            # Create list of terms that should not be reduced
            dont_reduce_terms = []
            for v in self.pos + self.neg:
                if not v in do_reduce_terms:
                    dont_reduce_terms.append(v)

            # Create expression that is not reduced
            not_reduced_expr = None
            if dont_reduce_terms and len(dont_reduce_terms) > 1:
                # Try to reduce the remaining terms that was not reduced at first
                not_reduced_expr = Sum(dont_reduce_terms).reduce_ops()
            elif dont_reduce_terms:
                not_reduced_expr = dont_reduce_terms[0]

            # Create return expression
            if not_reduced_expr:
                return Sum([reduced_expr, not_reduced_expr])
            return reduced_expr

        # Return self if we don't have any variables for which we can reduce
        # the sum
        return self

from floatvalue_obj import FloatValue
from symbol_obj     import Symbol
from product_obj    import Product
from fraction_obj   import Fraction

