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

class Sum(object):
#    __slots__ = ("val", "t", "pos", "neg", "expanded", "reduced", "_class", "_hash", "_repr", "_ops")
    __slots__ = ("val", "t", "pos", "neg", "expanded", "reduced", "_class", "_hash", "_repr")
#    __slots__ = ("val", "t", "pos", "neg", "reduced", "_class", "_hash", "_repr")
#    __slots__ = ("val", "t", "pos", "neg", "expanded", "reduced", "_class", "_hash")
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
#        pos_append = self.pos.append
#        neg_append = self.neg.append
        self.expanded = False
        self.reduced = False

        # Initialise class type
        self._class = "sum"

        if variables:
            # Remove nested Sums
            new_vars = []
#            new_vars_append = new_vars.append
            for v in variables:
#                if isinstance(v, Sum):
                if v._class == "sum":
                    new_vars += v.pos + v.neg
                    continue
                new_vars.append(v)
#                new_vars_append(v)

            floats = []
            # Loop variables and sort in positive and negative
            # Also collect all floats in 1 variable
            # We don't collect [x, x, x] into 3*x to avoid creating objects
            # Instead we do this when expanding the object
            for v in new_vars:
                # Skip zero terms
                if v.val ==  0.0:
                    continue
#                if isinstance(v, FloatValue):
                if v._class == "float":
                    floats.append(v)
                elif v.val < 0.0:
                    self.neg.append(v)
#                    neg_append(v)
                else:
                    self.pos.append(v)
#                    pos_append(v)

            # Only create new float if we have more than one. Also ignore it
            # if 1 + 1 - 2 = 0
            if len(floats) > 1:
                val = sum([v.val for v in floats])
                if val and val < 0.0:
#                    self.neg.append(FloatValue(val))
                    self.neg.append(create_float(val))
#                    neg_append(create_float(val))
                elif val:
#                    self.pos.append(FloatValue(val))
                    self.pos.append(create_float(val))
#                    pos_append(create_float(val))
            elif floats:
                if floats[0].val and floats[0].val < 0.0:
                    self.neg.append(floats[0])
#                    neg_append(floats[0])
                elif floats[0].val:
                    self.pos.append(floats[0])
#                    pos_append(floats[0])

        # If we don't have any variables the sum is zero
        else:
#            self.pos = [FloatValue(0)]
            self.pos = [create_float(0)]

        # Handle zero value
        if not self.pos + self.neg:
            self.val = 0.0
#            self.pos = [FloatValue(0)]
            self.pos = [create_float(0)]

        # Type is equal to the smallest type in both lists
        self.t = min([v.t for v in self.neg + self.pos])

        # Sort variables in positive and negative, (for representation)
        self.pos.sort()
        self.neg.sort()

        # Compute the representation now, such that we can use it directly
        # in the __eq__ and __ne__ methods (improves performance a bit, but
        # only when objects are cached).
#        self._repr = self.__repr__()
        self._repr = "Sum([%s])" % ", ".join([v._repr for v in self.pos + self.neg])

        # TODO: Use cache for hash instead
        self._hash = hash(self._repr)

    # Print functions
    def __repr__(self):
        "Representation for debugging"
#        return "Sum([%s])" % ", ".join([repr(v) for v in self.pos + self.neg])
#        if not self._repr:
#            self._repr = "Sum([%s])" % ", ".join([repr(v) for v in self.pos + self.neg])
        return self._repr

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

    # Hash (for lookup in {})
    def __hash__(self):
        "Use repr as hash"
#        if self._hash:
#            return self._hash
#        self._hash = hash(repr(self))
        return self._hash

    # Comparison
    def __eq__(self, other):
        "Two sums are equal if their list of variables are equal including sign"
#        if other and other._class in ("float", "sym", "prod", "sum", "frac"):
        if other:
            return self._repr == other._repr
#        if isinstance(other, Sum):
#        if other and other._class == "sum":
#            return self.pos == other.pos and self.neg == other.neg
        return False

    def __ne__(self, other):
        "Two sums are not equal if equal is false"
#        if other and other._class in ("float", "sym", "prod", "sum", "frac"):
        if other:
            return self._repr != other._repr
        return True

    def __lt__(self, other):
        # Symbols and products are always less
#        if isinstance(other, (FloatValue, Symbol, Product)):
        if other._class in ("float", "sym", "prod"):
            return False
        # Compare representation, to get members and sign correct
#        elif isinstance(other, Sum):
        elif other._class == "sum":
            return self.pos + self.neg < other.pos + other.neg
        # Fractions are always greater
        return True

    def __gt__(self, other):
        # Symbols and products are always less
#        if isinstance(other, (FloatValue, Symbol, Product)):
        if other._class in ("float", "sym", "prod"):
            return True
        # Compare representation, to get members and sign correct
#        elif isinstance(other, Sum):
        elif other._class == "sum":
            return self.pos + self.neg > other.pos + other.neg
        return False

    # Binary operators
    def __mul__(self, other):

        # If product will be zero
        if self.val == 0.0 or other.val == 0.0:
#            return FloatValue(0)
            return create_float(0)

        # NOTE: We expect expanded sub-expressions with no nested operators
        # Create list of new products using the '*' operator
        # TODO: Is this efficient?
        new_prods = [v*other for v in self.pos + self.neg]

        # Remove zero valued terms
        # TODO: Can this still happen?
        new_prods = [v for v in new_prods if v.val != 0.0]

        # Create new sum
        if not new_prods:
#            return FloatValue(0)
            return create_float(0)
        elif len(new_prods) > 1:
            # Expand sum to collect terms
#            return Sum(new_prods).expand()
            return create_sum(new_prods).expand()
        # TODO: Is it necessary to call expand?      
        return new_prods[0].expand()

    def __div__(self, other):

        # If division is illegal (this should definitely not happen)
        if other.val == 0.0:
            raise RuntimeError("Division by zero")

        # If fraction will be zero
        if self.val == 0.0:
#            return FloatValue(0)
            return create_float(0)

        # NOTE: assuming that we get expanded variables
        # If other is a Sum we can only return a fraction
        # TODO: We could check for equal sums if Sum.__eq__ could be trusted
        # As it is now (2*x + y) == (3*x + y), which works for the other things I do
        # NOTE: Expect that other is expanded i.e., x + x -> 2*x which can be handled
        # TODO: Fix (1 + y) / (x + x*y) -> 1 / x
        # Will this be handled when reducing operations on a fraction?
#        if isinstance(other, Sum):
        if other._class == "sum":
#            return Fraction(self, other)
            return create_fraction(self, other)


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
#            return FloatValue(0)
            return create_float(0)
        elif len(new_fracs) > 1:
#            return Sum(new_fracs)
            return create_sum(new_fracs)
        return new_fracs[0]

    # Public functions
    def ops(self):
        "Return number of operations to compute value of sum"
#        if self._ops:
#            return self._ops
        op = 0
        # Add the number of operations from sub-expressions
        for v in self.pos + self.neg:
            #  +1 for the +/- symbol
            op += v.ops() + 1

        # Subtract one operation as it only takes n-1 ops to sum n members
#        self._ops = op - 1
#        return self._ops
        return op - 1

    def expand(self):
        "Expand all members of the sum"

        if self.expanded:
            return self.expanded

        # TODO: This function might need some optimisation
        # Get all variables and expand (is OK because we don't have any nested Sums)
        variables = [v.expand() for v in self.pos + self.neg]

        # Remove nested sums that got created by expanding variables
        new_vars = []
#        new_append = new_vars.append
        for v in variables:
#            if isinstance(v, Sum):
            if v._class == "sum":
                new_vars += v.pos + v.neg
                continue
            new_vars.append(v)
#            new_append(v)

        # Sort variables into symbols, products and fractions (add floats
        # directly to new list, will be handled later). Add fractions if
        # possible else add to list.
        new_variables = []
#        append = new_variables.append
        syms = []
#        sym_append = syms.append
        prods = []
#        prod_append = prods.append
        frac_groups = {}
        # TODO: Rather than using '+', would it be more efficient to collect
        # the terms first?
        for v in new_vars:
#            print "v: ", repr(v)
#            if isinstance(v, FloatValue):
            # TODO: Should we also group fractions, or put this in a separate function?
#            if isinstance(v, (FloatValue, Fraction)):
            if v._class in ("float", "frac"):
                new_variables.append(v)
#                append(v)
#            elif isinstance(v, Symbol):
            elif v._class == "sym":
                syms.append(v)
#                sym_append(v)
#            elif isinstance(v, Product):
            elif v._class == "prod":
                prods.append(v)
#                prod_append(v)
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
#            if v._vrs in prod_groups:
                prod_groups[v.get_vrs()] += v
#                prod_groups[v._vrs] += v
            else:
                prod_groups[v.get_vrs()] = v
#                prod_groups[v._vrs] = v

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
        for k,v in sym_groups.iteritems():
            new_variables.append(v)
#            append(v)
        for k,v in prod_groups.iteritems():
            new_variables.append(v)
#            append(v)
#        for k,v in frac_groups.iteritems():
#            new_variables.append(v)
#            append(v)

        if len(new_variables) > 1:
            # Return new sum (will remove multiple instances of floats during construction)
#            return Sum(new_variables)
#            return create_sum(new_variables)
            self.expanded = create_sum(new_variables)
            return self.expanded
        elif new_variables:
            # If we just have one variable left, return it since it is already expanded
#            return new_variables[0]
            self.expanded = new_variables[0]
            return self.expanded
        raise RuntimeError("Where did the variables go?")

    def reduce_vartype(self, var_type):
        """Reduce expression with given var_type. It returns a list of tuples
        [(found, remain)], where 'found' is an expression that only has variables
        of type == var_type. If no variables are found, found=(). The 'remain'
        part contains the leftover after division by 'found' such that:
        self = Sum([f*r for f,r in self.reduce_vartype(Type)])."""

#        if var_type in self._reduce_vartype:
#            return self._reduce_vartype[var_type]

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
#        append = returns.append
        for f, r in found.iteritems():
            if len(r) > 1:
                # Use expand to group expressions
#                r = Sum(r).expand()
                r = create_sum(r).expand()
            elif r:
                r = r.pop()
            returns.append((f, r))
#            append((f, r))
        return returns
        self._reduce_vartype[var_type] = returns
        return returns

    def get_unique_vars(self, var_type):
        "Get unique variables (Symbols) as a set"
        # Loop all variables of self update the set
#        if var_type in self._unique_vars:
#            return self._unique_vars[var_type]
        var = set()
        update = var.update
        for v in self.pos + self.neg:
            update(v.get_unique_vars(var_type))
#        self._unique_vars[var_type] = var
        return var

    def get_var_occurrences(self):
        """Determine the number of minimum number of times all variables occurs
        in the expression. Returns a dictionary of variables and the number of
        times they occur. x*x + x returns {x:1}, x + y returns {}"""
        # NOTE: This function is only used if the numerator of a Fraction is a Sum

        # Get occurrences in first expression
        d0 = (self.pos + self.neg)[0].get_var_occurrences()
        for var in (self.pos + self.neg)[1:]:
            # Get the occurrences
            d = var.get_var_occurrences()
            # Delete those variables in d0 that are not in d
            for k, v in d0.items():
                if not k in d:
                    del d0[k]
            # Set the number of occurrences equal to the smallest number
            for k, v in d.iteritems():
                if k in d0:
                    d0[k] = min(d0[k], v)
        return d0

    def reduce_ops(self):
        "Reduce the number of operations needed to evaluate the sum"

        if self.reduced:
            return self.reduced
        # NOTE: Assuming that sum has already been expanded
        # TODO: Add test for this and handle case if it is not

        # TODO: The entire function looks expensive, can it be optimised?

        # TODO: It is not necessary to create a new Sum if we do not have more
        # than one Fraction
        # First group all fractions in the sum
        new_sum = group_fractions(self)
#        print "\nnew_sum: ", new_sum
#        if not isinstance(new_sum, Sum):
        if new_sum._class != "sum":
            self.reduced = new_sum.reduce_ops()
            return self.reduced
        # Loop all variables of the sum and collect the number of common
        # variables that can be factored out.
        common_vars = {}
        for var in new_sum.pos + new_sum.neg:

            # Get dictonary of occurrences
            # TODO: Don't save d, group with following line
            d = var.get_var_occurrences()
#            d = get_var_occurrences(var)

            # Add the variable and the number of occurrences to common dictionary
            for k, v in d.iteritems():
                if k in common_vars:
                    common_vars[k].append((v, var))
                else:
                    common_vars[k] = [(v, var)]

#        print "\ncommon_vars"
#        for k, v in common_vars.iteritems():
#            print "var: ", k
#            print "terms: ", Sum([t[1] for t in v])

        # Determine the maximum reduction for each variable
        # sorted as: {(x*x*y, x*y*z, 2*y):[2, [y]]}
        terms_reductions = {}
        for k, v in common_vars.iteritems():
            # If the number of expressions that can be reduced is only one
            # there is nothing to be done
            if len(v) > 1:
                # TODO: Is there a better way to compute the reduction gain
                # and the number of occurrences we should remove?

                # Get the list of number of occurences of 'k' in expressions
                # in 'v'
                occurrences = [t[0] for t in v]
                # Extract the terms of v
                terms = [t[1] for t in v]

                # Determine the favorable number of occurences and an estimate
                # of the maximum reduction for current variable
                fav_occur = 0
                reduc = 0
                for i in set(occurrences):
                    # Get number of terms that has a number of occcurences equal
                    # to or higher than the current number
                    num_terms = len([o for o in occurrences if o >= i])

                    # An estimate of the reduction in operations is:
                    # (number_of_terms - 1) * number_occurrences
                    new_reduc = (num_terms-1)*i
                    if new_reduc > reduc:
                        reduc = new_reduc
                        fav_occur = i

                # We need to reduce the expression with the favorable number of
                # occurrences of the current variable
                red_vars = [k]*fav_occur

                # If the list of terms is already present in the dictionary,
                # add the reduction count and the variables
                if tuple(terms) in terms_reductions:
                    terms_reductions[tuple(terms)][0] += reduc
                    terms_reductions[tuple(terms)][1] += red_vars
                else:
                    terms_reductions[tuple(terms)] = [reduc, red_vars]

        if terms_reductions:
#            print "\nterms_reductions"
#            for k, v in terms_reductions.iteritems():
#                print Sum(k)
#                print v[0], v[1]

            # Invert dictionary of terms
            reductions_terms = dict([((v[0], tuple(v[1])), k) for k, v in terms_reductions.iteritems()])

            # Create a sorted list of those variables that give the highest
            # reduction
            sorted_reduc_var = [k for k, v in reductions_terms.iteritems()]
            sorted_reduc_var.sort(lambda x, y: cmp(x[0], y[0]))
            sorted_reduc_var.reverse()

            # Create a new dictionary of terms that should be reduced, if some
            # terms overlap, only pick the one which give the highest reduction to
            # ensure that a*x*x + b*x*x + x*x*y + 2*y -> x*x*(a + b + y) + 2*y NOT 
            # x*x*(a + b) + y*(2 + x*x)
            reduction_vars = {}
            rejections = {}
            for var in sorted_reduc_var:
                terms = reductions_terms[var]
                if overlap(terms, reduction_vars) or overlap(terms, rejections):
                    rejections[var[1]] = terms
                else:
                    reduction_vars[var[1]] = terms

#            print "\nreduction_vars"
#            for k, v in reduction_vars.iteritems():
#                print k
#                print Sum(v)

            # Reduce each set of terms with appropriate variables
            all_reduced_terms = []
            reduced_expressions = []
            for reduc_var, terms in reduction_vars.iteritems():
#                print "reduc_var: ", reduc_var
#                print "terms: ", Sum(terms)
                # Add current terms to list of all variables that have been reduced
                all_reduced_terms += list(terms)

                # Create variable that we will use to reduce the terms
                reduction_var = None
                if len(reduc_var) > 1:
#                    reduction_var = Product(reduc_var)
                    reduction_var = create_product(list(reduc_var))
                else:
                    reduction_var = reduc_var[0]

                # Reduce all terms that need to be reduced
                reduced_terms = [t.reduce_var(reduction_var) for t in terms]

                # Create reduced expression
                reduced_expr = None
                if len(reduced_terms) > 1:
                    # Try to reduce the reduced terms further
#                    print "recursion on reduced terms"
#                    reduced_expr = Product([reduction_var, Sum(reduced_terms).reduce_ops()])
#                    reduced_expr = create_product([reduction_var, Sum(reduced_terms).reduce_ops()])
                    reduced_expr = create_product([reduction_var, create_sum(reduced_terms).reduce_ops()])
                else:
#                    reduced_expr = Product(reduction_var, reduced_terms[0])
                    reduced_expr = create_product(reduction_var, reduced_terms[0])

                # Add reduced expression to list of reduced expressions
                reduced_expressions.append(reduced_expr)

            # Create list of terms that should not be reduced
            dont_reduce_terms = []
#            dont_reduce_terms_append = dont_reduce_terms.append
            for v in new_sum.pos + new_sum.neg:
                if not v in all_reduced_terms:
                    dont_reduce_terms.append(v)
#                    dont_reduce_terms_append(v)
#            print "dont reduce: ", dont_reduce_terms
#            print "reduced expr: ", reduced_expressions

            # Create expression from terms that was not reduced
            not_reduced_expr = None
            if dont_reduce_terms and len(dont_reduce_terms) > 1:
#                print "recursion on not reduced terms"
                # Try to reduce the remaining terms that were not reduced at first
#                print "expr: ", Sum(dont_reduce_terms)
#                not_reduced_expr = Sum(dont_reduce_terms).reduce_ops()
                not_reduced_expr = create_sum(dont_reduce_terms).reduce_ops()
            elif dont_reduce_terms:
                not_reduced_expr = dont_reduce_terms[0]

            # Create return expression
            if not_reduced_expr:
#                return Sum(reduced_expressions + [not_reduced_expr])
#                return create_sum(reduced_expressions + [not_reduced_expr])
                self.reduced = create_sum(reduced_expressions + [not_reduced_expr])
            elif len(reduced_expressions) > 1:
#                return Sum(reduced_expressions)
#                return create_sum(reduced_expressions)
                self.reduced = create_sum(reduced_expressions)
            else:
                self.reduced = reduced_expressions[0]
            return self.reduced

        # Return self if we don't have any variables for which we can reduce
        # the sum
        self.reduced = self
        return self.reduced

def overlap(l, d):
    "Check if a member in list l is in the value (list) of dictionary d"
    for m in l:
        for k, v in d.iteritems():
            if m in v:
                return True
    return False

def group_fractions(expr):
    "Group Fractions in a Sum: 2/x + y/x -> (2 + y)/x"
#    if not isinstance(expr, Sum):
    if expr._class != "sum":
        return expr

    # Loop variables and group those with common denominator
    not_frac = []
#    not_frac_append = not_frac.append
    fracs = {}
    for v in expr.pos + expr.neg:
#        if isinstance(v, Fraction):
        if v._class == "frac":
            if v.denom in fracs:
                fracs[v.denom][1].append(v.num)
                fracs[v.denom][0] += 1
            else:
                fracs[v.denom] = [1, [v.num], v]
            continue
        not_frac.append(v)
#        not_frac_append(v)
    if not fracs:
        return expr

    for k, v in fracs.iteritems():
        if v[0] > 1:
            # TODO: Is it possible to avoid expanding the Sum?
            # I think we have to because x/a + 2*x/a -> 3*x/a
#            not_frac.append(Fraction(Sum(v[1]).expand(), k))
#            not_frac.append(Fraction(create_sum(v[1]).expand(), k))
            not_frac.append(create_fraction(create_sum(v[1]).expand(), k))
#            not_frac_append(create_fraction(create_sum(v[1]).expand(), k))
        else:
            not_frac.append(v[2])
#            not_frac_append(v[2])

    if len(not_frac) > 1:
#        return Sum(not_frac)
        return create_sum(not_frac)
    return not_frac[0]


from floatvalue_obj import FloatValue
from symbol_obj     import Symbol
from product_obj    import Product
from fraction_obj   import Fraction

