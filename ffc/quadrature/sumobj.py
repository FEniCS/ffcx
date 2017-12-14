# -*- coding: utf-8 -*-
"This file implements a class to represent a sum."

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

from ufl.utils.sorting import sorted_by_key

# FFC modules.
from ffc.log import error
from ffc.quadrature.cpp import format

# FFC quadrature modules.
from .symbolics import create_float
from .symbolics import create_product
from .symbolics import create_sum
from .symbolics import create_fraction
from .expr import Expr

from .floatvalue import FloatValue


class Sum(Expr):
    __slots__ = ("vrs", "_expanded", "_reduced")

    def __init__(self, variables):
        """Initialise a Sum object, it derives from Expr and contains the
        additional variables:

        vrs       - list, a list of variables.
        _expanded - object, an expanded object of self, e.g.,
                    self = 'x + x'-> self._expanded = 2*x (a product).
        _reduced  - object, a reduced object of self, e.g.,
                    self = '2*x + x*y'-> self._reduced = x*(2 + y) (a product).
        NOTE: self._prec = 3."""

        # Initialise value, list of variables, class, expanded and reduced.
        self.val = 1.0
        self.vrs = []
        self._prec = 3
        self._expanded = False
        self._reduced = False

        # Get epsilon
        EPS = format["epsilon"]
        # Process variables if we have any.
        if variables:
            # Loop variables and remove nested Sums and collect all floats in
            # 1 variable. We don't collect [x, x, x] into 3*x to avoid creating
            # objects, instead we do this when expanding the object.
            float_val = 0.0
            for var in variables:
                # Skip zero terms.
                if abs(var.val) < EPS:
                    continue
                elif var._prec == 0:  # float
                    float_val += var.val
                    continue
                elif var._prec == 3:  # sum
                    # Loop and handle variables of nested sum.
                    for v in var.vrs:
                        if abs(v.val) < EPS:
                            continue
                        elif v._prec == 0:  # float
                            float_val += v.val
                            continue
                        self.vrs.append(v)
                    continue
                self.vrs.append(var)

            # Only create new float if value is different from 0.
            if abs(float_val) > EPS:
                self.vrs.append(create_float(float_val))

        # If we don't have any variables the sum is zero.
        else:
            self.val = 0.0
            self.vrs = [create_float(0)]

        # Handle zero value.
        if not self.vrs:
            self.val = 0.0
            self.vrs = [create_float(0)]

        # Type is equal to the smallest type in both lists.
        self.t = min([v.t for v in self.vrs])

        # Sort variables, (for representation).
        self.vrs.sort()

        # Compute the representation now, such that we can use it directly
        # in the __eq__ and __ne__ methods (improves performance a bit, but
        # only when objects are cached).
        self._repr = "Sum([%s])" % ", ".join([v._repr for v in self.vrs])

        # Use repr as hash value.
        self._hash = hash(self._repr)

    # Print functions.
    def __str__(self):
        "Simple string representation which will appear in the generated code."
        # First add all the positive variables using plus, then add all
        # negative variables.
        s = format["add"]([str(v) for v in self.vrs if not v.val < 0]) +\
            "".join([str(v) for v in self.vrs if v.val < 0])
        # Group only if we have more that one variable.
        if len(self.vrs) > 1:
            return format["grouping"](s)
        return s

    # Binary operators.
    def __add__(self, other):
        "Addition by other objects."
        # Return a new sum
        return create_sum([self, other])

    def __sub__(self, other):
        "Subtract other objects."
        # Return a new sum
        return create_sum([self, create_product([FloatValue(-1), other])])

    def __mul__(self, other):
        "Multiplication by other objects."
        # If product will be zero.
        if self.val == 0.0 or other.val == 0.0:
            return create_float(0)

        # NOTE: We expect expanded sub-expressions with no nested operators.
        # Create list of new products using the '*' operator
        # TODO: Is this efficient?
        new_prods = [v * other for v in self.vrs]

        # Remove zero valued terms.
        # TODO: Can this still happen?
        new_prods = [v for v in new_prods if v.val != 0.0]

        # Create new sum.
        if not new_prods:
            return create_float(0)
        elif len(new_prods) > 1:
            # Expand sum to collect terms.
            return create_sum(new_prods).expand()
        # TODO: Is it necessary to call expand?
        return new_prods[0].expand()

    def __truediv__(self, other):
        "Division by other objects."
        # If division is illegal (this should definitely not happen).
        if other.val == 0.0:
            error("Division by zero.")

        # If fraction will be zero.
        if self.val == 0.0:
            return create_float(0)

        # NOTE: assuming that we get expanded variables.
        # If other is a Sum we can only return a fraction.
        # TODO: We could check for equal sums if Sum.__eq__ could be trusted.
        # As it is now (2*x + y) == (3*x + y), which works for the other things I do.
        # NOTE: Expect that other is expanded i.e., x + x -> 2*x which can be handled.
        # TODO: Fix (1 + y) / (x + x*y) -> 1 / x
        # Will this be handled when reducing operations on a fraction?
        if other._prec == 3:  # sum
            return create_fraction(self, other)

        # NOTE: We expect expanded sub-expressions with no nested operators.
        # Create list of new products using the '*' operator.
        # TODO: Is this efficient?
        new_fracs = [v / other for v in self.vrs]

        # Remove zero valued terms.
        # TODO: Can this still happen?
        new_fracs = [v for v in new_fracs if v.val != 0.0]

        # Create new sum.
        # TODO: No need to call expand here, using the '/' operator should have
        # taken care of this.
        if not new_fracs:
            return create_float(0)
        elif len(new_fracs) > 1:
            return create_sum(new_fracs)
        return new_fracs[0]

    __div__ = __truediv__

    # Public functions.
    def expand(self):
        "Expand all members of the sum."

        # If sum is already expanded, simply return the expansion.
        if self._expanded:
            return self._expanded

        # TODO: This function might need some optimisation.

        # Sort variables into symbols, products and fractions (add
        # floats directly to new list, will be handled later). Add
        # fractions if possible else add to list.
        new_variables = []
        syms = []
        prods = []
        # TODO: Rather than using '+', would it be more efficient to
        # collect the terms first?
        for var in self.vrs:
            exp = var.expand()
            # TODO: Should we also group fractions, or put this in a
            # separate function?
            if exp._prec in (0, 4):  # float or frac
                new_variables.append(exp)
            elif exp._prec == 1:  # sym
                syms.append(exp)
            elif exp._prec == 2:  # prod
                prods.append(exp)
            elif exp._prec == 3:  # sum
                for v in exp.vrs:
                    if v._prec in (0, 4):  # float or frac
                        new_variables.append(v)
                    elif v._prec == 1:  # sym
                        syms.append(v)
                    elif v._prec == 2:  # prod
                        prods.append(v)

        # Sort all variables in groups: [2*x, -7*x], [(x + y), (2*x +
        # 4*y)] etc.  First handle product in order to add symbols if
        # possible.
        prod_groups = {}
        for v in prods:
            if v.get_vrs() in prod_groups:
                prod_groups[v.get_vrs()] += v
            else:
                prod_groups[v.get_vrs()] = v

        sym_groups = {}
        # Loop symbols and add to appropriate groups.
        for v in syms:
            # First try to add to a product group.
            if (v,) in prod_groups:
                prod_groups[(v,)] += v
            # Then to other symbols.
            elif v in sym_groups:
                sym_groups[v] += v
            # Create a new entry in the symbols group.
            else:
                sym_groups[v] = v

        # Loop groups and add to new variable list.
        for k, v in sorted_by_key(sym_groups):
            new_variables.append(v)
        for k, v in sorted_by_key(prod_groups):
            new_variables.append(v)

        if len(new_variables) > 1:
            # Return new sum (will remove multiple instances of floats
            # during construction).
            self._expanded = create_sum(sorted(new_variables))
            return self._expanded
        elif new_variables:
            # If we just have one variable left, return it since it is
            # already expanded.
            self._expanded = new_variables[0]
            return self._expanded
        error("Where did the variables go?")

    def get_unique_vars(self, var_type):
        "Get unique variables (Symbols) as a set."
        # Loop all variables of self update the set.
        var = set()
        for v in self.vrs:
            var.update(v.get_unique_vars(var_type))
        return var

    def get_var_occurrences(self):
        """Determine the number of minimum number of times all variables
        occurs in the expression. Returns a dictionary of variables
        and the number of times they occur. x*x + x returns {x:1}, x +
        y returns {}.

        """
        # NOTE: This function is only used if the numerator of a
        # Fraction is a Sum.

        # Get occurrences in first expression.
        d0 = self.vrs[0].get_var_occurrences()
        for var in self.vrs[1:]:
            # Get the occurrences.
            d = var.get_var_occurrences()
            # Delete those variables in d0 that are not in d.
            for k, v in list(d0.items()):
                if k not in d:
                    del d0[k]
            # Set the number of occurrences equal to the smallest
            # number.
            for k, v in sorted_by_key(d):
                if k in d0:
                    d0[k] = min(d0[k], v)
        return d0

    def ops(self):
        "Return number of operations to compute value of sum."
        # Subtract one operation as it only takes n-1 ops to sum n
        # members.
        op = -1

        # Add the number of operations from sub-expressions.
        for v in self.vrs:
            #  +1 for the +/- symbol.
            op += v.ops() + 1
        return op

    def reduce_ops(self):
        "Reduce the number of operations needed to evaluate the sum."

        if self._reduced:
            return self._reduced
        # NOTE: Assuming that sum has already been expanded.
        # TODO: Add test for this and handle case if it is not.

        # TODO: The entire function looks expensive, can it be optimised?

        # TODO: It is not necessary to create a new Sum if we do not
        # have more than one Fraction.

        # First group all fractions in the sum.
        new_sum = _group_fractions(self)
        if new_sum._prec != 3:  # sum
            self._reduced = new_sum.reduce_ops()
            return self._reduced
        # Loop all variables of the sum and collect the number of
        # common variables that can be factored out.
        common_vars = {}
        for var in new_sum.vrs:
            # Get dictonary of occurrences and add the variable and
            # the number of occurrences to common dictionary.
            for k, v in sorted_by_key(var.get_var_occurrences()):
                if k in common_vars:
                    common_vars[k].append((v, var))
                else:
                    common_vars[k] = [(v, var)]

        # Determine the maximum reduction for each variable sorted as:
        # {(x*x*y, x*y*z, 2*y):[2, [y]]}.
        terms_reductions = {}
        for k, v in sorted_by_key(common_vars):
            # If the number of expressions that can be reduced is only
            # one there is nothing to be done.
            if len(v) > 1:
                # TODO: Is there a better way to compute the reduction
                # gain and the number of occurrences we should remove?

                # Get the list of number of occurences of 'k' in
                # expressions in 'v'.
                occurrences = [t[0] for t in v]

                # Determine the favorable number of occurences and an
                # estimate of the maximum reduction for current
                # variable.
                fav_occur = 0
                reduc = 0
                for i in set(occurrences):
                    # Get number of terms that has a number of
                    # occcurences equal to or higher than the current
                    # number.
                    num_terms = len([o for o in occurrences if o >= i])

                    # An estimate of the reduction in operations is:
                    # (number_of_terms - 1) * number_occurrences.
                    new_reduc = (num_terms - 1) * i
                    if new_reduc > reduc:
                        reduc = new_reduc
                        fav_occur = i

                # Extract the terms of v where the number of
                # occurrences is equal to or higher than the most
                # favorable number of occurrences.
                terms = sorted([t[1] for t in v if t[0] >= fav_occur])

                # We need to reduce the expression with the favorable
                # number of occurrences of the current variable.
                red_vars = [k] * fav_occur

                # If the list of terms is already present in the
                # dictionary, add the reduction count and the
                # variables.
                if tuple(terms) in terms_reductions:
                    terms_reductions[tuple(terms)][0] += reduc
                    terms_reductions[tuple(terms)][1] += red_vars
                else:
                    terms_reductions[tuple(terms)] = [reduc, red_vars]

        if terms_reductions:
            # Invert dictionary of terms.
            reductions_terms = dict([((v[0], tuple(v[1])), k) for k,
                                     v in terms_reductions.items()])

            # Create a sorted list of those variables that give the
            # highest reduction.
            sorted_reduc_var = sorted(reductions_terms.keys(),
                                      reverse=True)

            # Create a new dictionary of terms that should be reduced,
            # if some terms overlap, only pick the one which give the
            # highest reduction to ensure that a*x*x + b*x*x + x*x*y +
            # 2*y -> x*x*(a + b + y) + 2*y NOT x*x*(a + b) + y*(2 +
            # x*x).
            reduction_vars = {}
            rejections = {}
            for var in sorted_reduc_var:
                terms = reductions_terms[var]
                if _overlap(terms, reduction_vars) or _overlap(terms,
                                                               rejections):
                    rejections[var[1]] = terms
                else:
                    reduction_vars[var[1]] = terms

            # Reduce each set of terms with appropriate variables.
            all_reduced_terms = []
            reduced_expressions = []
            for reduc_var, terms in sorted(reduction_vars.items()):

                # Add current terms to list of all variables that have
                # been reduced.
                all_reduced_terms += list(terms)

                # Create variable that we will use to reduce the terms.
                reduction_var = None
                if len(reduc_var) > 1:
                    reduction_var = create_product(list(reduc_var))
                else:
                    reduction_var = reduc_var[0]

                # Reduce all terms that need to be reduced.
                reduced_terms = [t.reduce_var(reduction_var) for t in terms]

                # Create reduced expression.
                reduced_expr = None
                if len(reduced_terms) > 1:
                    # Try to reduce the reduced terms further.
                    reduced_expr = create_product([reduction_var,
                                                   create_sum(reduced_terms).reduce_ops()])
                else:
                    reduced_expr = create_product(reduction_var,
                                                  reduced_terms[0])

                # Add reduced expression to list of reduced
                # expressions.
                reduced_expressions.append(reduced_expr)

            # Create list of terms that should not be reduced.
            dont_reduce_terms = []
            for v in new_sum.vrs:
                if v not in all_reduced_terms:
                    dont_reduce_terms.append(v)

            # Create expression from terms that was not reduced.
            not_reduced_expr = None
            if dont_reduce_terms and len(dont_reduce_terms) > 1:
                # Try to reduce the remaining terms that were not
                # reduced at first.
                not_reduced_expr = create_sum(dont_reduce_terms).reduce_ops()
            elif dont_reduce_terms:
                not_reduced_expr = dont_reduce_terms[0]

            # Create return expression.
            if not_reduced_expr:
                self._reduced = create_sum(reduced_expressions + [not_reduced_expr])
            elif len(reduced_expressions) > 1:
                self._reduced = create_sum(reduced_expressions)
            else:
                self._reduced = reduced_expressions[0]

            return self._reduced

        # Return self if we don't have any variables for which we can
        # reduce the sum.
        self._reduced = self
        return self._reduced

    def reduce_vartype(self, var_type):
        """Reduce expression with given var_type. It returns a list of tuples
        [(found, remain)], where 'found' is an expression that only
        has variables of type == var_type. If no variables are found,
        found=(). The 'remain' part contains the leftover after
        division by 'found' such that: self = Sum([f*r for f,r in
        self.reduce_vartype(Type)]).

        """
        found = {}
        # Loop members and reduce them by vartype.
        for v in self.vrs:
            for f, r in v.reduce_vartype(var_type):
                if f in found:
                    found[f].append(r)
                else:
                    found[f] = [r]

        # Create the return value.
        returns = []
        for f, r in sorted_by_key(found):
            if len(r) > 1:
                # Use expand to group expressions.
                r = create_sum(r)
            elif r:
                r = r.pop()
            returns.append((f, r))
        return sorted(returns)


def _overlap(l, d):
    "Check if a member in list l is in the value (list) of dictionary d."
    for m in l:
        for k, v in sorted_by_key(d):
            if m in v:
                return True
    return False


def _group_fractions(expr):
    "Group Fractions in a Sum: 2/x + y/x -> (2 + y)/x."
    if expr._prec != 3:  # sum
        return expr

    # Loop variables and group those with common denominator.
    not_frac = []
    fracs = {}
    for v in expr.vrs:
        if v._prec == 4:  # frac
            if v.denom in fracs:
                fracs[v.denom][1].append(v.num)
                fracs[v.denom][0] += 1
            else:
                fracs[v.denom] = [1, [v.num], v]
            continue
        not_frac.append(v)
    if not fracs:
        return expr

    # Loop all fractions and create new ones using an appropriate
    # numerator.
    for k, v in sorted(fracs.items()):
        if v[0] > 1:
            # TODO: Is it possible to avoid expanding the Sum?  I
            # think we have to because x/a + 2*x/a -> 3*x/a.
            not_frac.append(create_fraction(create_sum(v[1]).expand(), k))
        else:
            not_frac.append(v[2])

    # Create return value.
    if len(not_frac) > 1:
        return create_sum(not_frac)
    return not_frac[0]
