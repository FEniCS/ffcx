# -*- coding: utf-8 -*-
"Some simple functions for manipulating expressions symbolically"

# Copyright (C) 2008-2010 Kristian B. Oelgaard
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

# FFC modules
from ffc.log import error

from collections import deque


def split_expression(expression, format, operator, allow_split=False):
    """Split the expression at the given operator, return list.
    Do not split () or [] unless told to split (). This is to enable easy count
    of double operations which can be in (), but in [] we only have integer operations."""

    # Get formats
    access = format["component"]("", [""])
    group = format["grouping"]("")
    la = access[0]
    ra = access[1]
    lg = group[0]
    rg = group[1]

    # Split with given operator
    prods = deque(expression.split(operator))
    new_prods = [prods.popleft()]

    while prods:
        # Continue while we still have list of potential products p is
        # the first string in the product
        p = prods.popleft()
        # If the number of "[" and "]" doesn't add up in the last
        # entry of the new_prods list, add p and see if it helps for
        # next iteration
        if new_prods[-1].count(la) != new_prods[-1].count(ra):
            new_prods[-1] = operator.join([new_prods[-1], p])
        # If the number of "(" and ")" doesn't add up (and we didn't
        # allow a split) in the last entry of the new_prods list, add
        # p and see if it helps for next iteration
        elif new_prods[-1].count(lg) != new_prods[-1].count(rg) and not allow_split:
            new_prods[-1] = operator.join([new_prods[-1], p])
        # If everything was fine, we can start a new entry in the
        # new_prods list
        else:
            new_prods.append(p)

    return new_prods


def operation_count(expression, format):
    """This function returns the number of double operations in an
    expression.  We do split () but not [] as we only have unsigned
    integer operations in [].

    """

    # Note we do not subtract 1 for the additions, because there is
    # also an assignment involved
    adds = len(split_expression(expression, format, format["add"](["", ""]),
                                True)) - 1
    mults = len(split_expression(expression, format,
                                 format["multiply"](["", ""]), True)) - 1
    return mults + adds


def get_simple_variables(expression, format):
    """This function takes as argument an expression (preferably expanded):
      expression = "x*x + y*x + x*y*z"

    returns a list of products and a dictionary:

      prods = ["x*x", "y*x", "x*y*z"]
      variables = {variable: [num_occurences, [pos_in_prods]]}
      variables = {"x":[3, [0,1,2]], "y":[2, [1,2]], "z":[1, [2]]}

    """

    # Get formats
    add = format["add"](["", ""])
    mult = format["multiply"](["", ""])
    format_float = format["floating point"]

    prods = split_expression(expression, format, add)
    prods = [p for p in prods if p]

    variables = {}
    for i, p in enumerate(prods):
        # Only extract unique variables
        vrs = list(set(split_expression(p, format, mult)))
        for v in vrs:
            # Try to convert variable to floats and back (so '2' == '2.0' etc.)
            try:
                v = format_float(float(v))
            except Exception:
                pass
            if v in variables:
                variables[v][0] += 1
                variables[v][1].append(i)
            else:
                variables[v] = [1, [i]]
    return (prods, variables)


def group_vars(expr, format):
    """Group variables in an expression, such that:
    "x + y + z + 2*y + 6*z" = "x + 3*y + 7*z"
    "x*x + x*x + 2*x + 3*x + 5" = "2.0*x*x + 5.0*x + 5"
    "x*y + y*x + 2*x*y + 3*x + 0*x + 5" = "5.0*x*y + 3.0*x + 5"
    "(y + z)*x + 5*(y + z)*x" = "6.0*(y + z)*x"
    "1/(x*x) + 2*1/(x*x) + std::sqrt(x) + 6*std::sqrt(x)" = "3*1/(x*x) + 7*std::sqrt(x)"
    """

    # Get formats
    format_float = format["floating point"]
    add = format["add"](["", ""])
    mult = format["multiply"](["", ""])

    new_prods = {}

    # Get list of products
    prods = split_expression(expr, format, add)

    # Loop products and collect factors
    for p in prods:
        # Get list of variables, and do a basic sort
        vrs = split_expression(p, format, mult)
        factor = 1
        new_var = []

        # Try to multiply factor with variable, else variable must be
        # multiplied by factor later
        # If we don't have a variable, set factor to zero and break
        for v in vrs:
            if v:
                try:
                    f = float(v)
                    factor *= f
                except Exception:
                    new_var.append(v)
            else:
                factor = 0
                break

        # Create new variable that must be multiplied with factor. Add
        # this variable to dictionary, if it already exists add factor
        # to other factors
        new_var.sort()
        new_var = mult.join(new_var)
        if new_var in new_prods:
            new_prods[new_var] += factor
        else:
            new_prods[new_var] = factor

    # Reset products
    prods = []
    for prod, f in sorted_by_key(new_prods):
        # If we have a product append mult of both
        if prod:
            # If factor is 1.0 we don't need it
            if f == 1.0:
                prods.append(prod)
            else:
                prods.append(mult.join([format_float(f), prod]))
        # If we just have a factor
        elif f:
            prods.append(format_float(f))

    prods.sort()
    return add.join(prods)


def reduction_possible(variables):
    """Find the variable that occurs in the most products, if more
    variables occur the same number of times and in the same products
    add them to list.

    """

    # Find the variable that appears in the most products
    max_val = 1
    max_var = ""
    max_vars = []
    for key, val in sorted_by_key(variables):
        if max_val < val[0]:
            max_val = val[0]
            max_var = key

    # If we found a variable that appears in products multiple times,
    # check if other variables appear in the exact same products
    if max_var:
        for key, val in sorted_by_key(variables):
            # Check if we have more variables in the same products
            if max_val == val[0] and variables[max_var][1] == val[1]:
                max_vars.append(key)
    return max_vars


def is_constant(variable, format, constants=[], from_is_constant=False):
    """Determine if a variable is constant or not.  The function accepts
    an optional list of variables (loop indices) that will be regarded
    as constants for the given variable. If none are supplied it is
    assumed that all array accesses will result in a non-constant
    variable.

    v = 2.0,          is constant
    v = Jinv_00*det,  is constant
    v = w[0][1],      is constant
    v = 2*w[0][1],    is constant
    v = W0[ip],       is constant if constants = ['ip'] else not
    v = P_t0[ip][j],  is constant if constants = ['j','ip'] else not

    """

    # Get formats
    access = format["array access"]("")
    add = format["add"](["", ""])
    mult = format["multiply"](["", ""])

    l = access[0]  # noqa: E741
    r = access[1]

    if not variable.count(l) == variable.count(r):
        print("variable: ", variable)
        error("Something wrong with variable")

    # Be sure that we don't have a compound
    variable = expand_operations(variable, format)

    prods = split_expression(variable, format, add)

    # Loop all products and variables and check if they're constant
    for p in prods:
        vrs = split_expression(p, format, mult)
        for v in vrs:
            # Check if each variable is constant, if just one fails
            # the entire variable is considered not to be constant
            const_var = False

            # If variable is in constants, well....
            if v in constants:
                const_var = True
                continue

            # If we don't have any '[' or ']' we have a constant
            # (unless we're dealing with a call from this funtions)
            elif not v.count(l) and not from_is_constant:
                const_var = True
                continue

            # If we have an array access variable, see if the index is
            # regarded a constant
            elif v.count(l):

                # Check if access is OK ('[' is before ']')
                if not v.index(l) < v.index(r):
                    print("variable: ", v)
                    error("Something is wrong with the array access")

                # Auxiliary variables
                index = ""
                left = 0
                inside = False
                indices = []

                # Loop all characters in variable and find indices
                for c in v:

                    # If character is ']' reduce left count
                    if c == r:
                        left -= 1

                    # If the '[' count has returned to zero, we have a
                    # complete index
                    if left == 0 and inside:
                        const_index = False  # Aux. var
                        if index in constants:
                            const_index = True

                        try:
                            int(index)
                            const_index = True
                        except Exception:
                            # Last resort, call recursively
                            if is_constant(index, format, constants, True):
                                const_index = True
                            pass

                        # Append index and reset values
                        if const_index:
                            indices.append(const_index)
                        else:
                            indices = [False]
                            break
                        index = ""
                        inside = False

                    # If we're inside an access, add character to index
                    if inside:
                        index += c

                    # If character is '[' increase the count, and
                    # we're inside an access
                    if c == l:
                        inside = True
                        left += 1

                # If all indices were constant, the variable is constant
                if all(indices):
                    const_var = True
                    continue

            else:
                # If it is a float, it is also constant
                try:
                    float(v)
                    const_var = True
                    continue
                except Exception:
                    pass

            # I no tests resulted in a constant variable, there is no
            # need to continue
            if not const_var:
                return False

    # If all variables were constant return True
    return True


def expand_operations(expression, format):
    """This function expands an expression and returns the value. E.g.,
    ((x + y))             --> x + y
    2*(x + y)             --> 2*x + 2*y
    (x + y)*(x + y)       --> x*x + y*y + 2*x*y
    z*(x*(y + 3) + 2) + 1 --> 1 + 2*z + x*y*z + x*z*3
    z*((y + 3)*x + 2) + 1 --> 1 + 2*z + x*y*z + x*z*3"""

    # Get formats
    add = format["add"](["", ""])
    mult = format["multiply"](["", ""])
    group = format["grouping"]("")
    l = group[0]  # noqa: E741
    r = group[1]

    # Check that we have the same number of left/right parenthesis in
    # expression
    if not expression.count(l) == expression.count(r):
        error("Number of left/right parenthesis do not match")

    # If we don't have any parenthesis, group variables and return
    if expression.count(l) == 0:
        return group_vars(expression, format)

    # Get list of additions
    adds = split_expression(expression, format, add)
    new_adds = []

    # Loop additions and get products
    for a in adds:
        prods = sorted(split_expression(a, format, mult))
        new_prods = []

        # FIXME: Should we use deque here?
        expanded = []
        for i, p in enumerate(prods):
            # If we have a group, expand inner expression
            if p[0] == l and p[-1] == r:
                # Add remaining products to new products and multiply
                # with all terms from expanded variable
                expanded_var = expand_operations(p[1:-1], format)
                expanded.append(split_expression(expanded_var, format, add))

            # Else, just add variable to list of new products
            else:
                new_prods.append(p)

        if expanded:
            # Combine all expanded variables and multiply by factor
            while len(expanded) > 1:
                first = expanded.pop(0)
                second = expanded.pop(0)
                expanded = [[mult.join([i] + [j]) for i in first for j in second]] + expanded
            new_adds += [mult.join(new_prods + [e]) for e in expanded[0]]
        else:
            # Else, just multiply products and add to list of products
            new_adds.append(mult.join(new_prods))

    # Group variables and return
    return group_vars(add.join(new_adds), format)


def reduce_operations(expression, format):
    """This function reduces the number of opertions needed to compute a
    given expression. It looks for the variable that appears the most
    and groups terms containing this variable inside parenthesis. The
    function is called recursively until no further reductions are
    possible.

    "x + y + x" = 2*x + y
    "x*x + 2.0*x*y + y*y" = y*y + (2.0*y + x)*x, not (x + y)*(x + y) as it should be!!
    z*x*y + z*x*3 + 2*z + 1" = z*(x*(y + 3) + 2) + 1

    """

    # Get formats
    add = format["add"](["", ""])
    mult = format["multiply"](["", ""])

    # Be sure that we have an expanded expression
    expression = expand_operations(expression, format)

    # Group variables to possibly reduce complexity
    expression = group_vars(expression, format)

    # Get variables and products
    prods, variables = get_simple_variables(expression, format)

    # Get the variables for which we can reduce the expression
    max_vars = reduction_possible(variables)
    new_prods = []
    no_mult = []
    max_vars.sort()

    # If we have variables that can be moved outside
    if max_vars:
        for p in prods:
            # Get the list of variables in current product
            li = sorted(split_expression(p, format, mult))

            # If the list of products is the same as what we intend of
            # moving outside the parenthesis, leave it (because x +
            # x*x + x*y should be x + (x + y)*x NOT (1.0 + x + y)*x)
            if li == max_vars:
                no_mult.append(p)
                continue
            else:
                # Get list of all variables from max_vars that are in
                # li
                indices = [i for i in max_vars if i in li]
                # If not all were present add to list of terms that
                # shouldn't be multiplied with variables and continue
                if indices != max_vars:
                    no_mult.append(p)
                    continue

            # Remove variables that we are moving outside
            for v in max_vars:
                li.remove(v)

            # Add to list of products
            p = mult.join(li)
            new_prods.append(p)

        # Sort lists
        no_mult.sort()
        new_prods.sort()
    else:
        # No reduction possible
        return expression

    # Recursively reduce sums with and without reduced variable
    new_prods = add.join(new_prods)
    if new_prods:
        new_prods = reduce_operations(new_prods, format)
    if no_mult:
        no_mult = [reduce_operations(add.join(no_mult), format)]

    # Group new products if we have a sum
    g = new_prods
    len_new_prods = len(split_expression(new_prods, format, add))
    if len_new_prods > 1:
        g = format["grouping"](new_prods)

    # The new expression is the sum of terms that couldn't be reduced
    # and terms that could be reduced multiplied by the reduction
    # e.g., expr = z + (x + y)*x
    new_expression = add.join(no_mult + [mult.join([g, mult.join(max_vars)])])

    return new_expression


def get_geo_terms(expression, geo_terms, offset, format):
    """This function returns a new expression where all geometry terms
    have been substituted with geometry declarations, these
    declarations are added to the geo_terms dictionary.

    """

    # Get formats
    add = format["add"](["", ""])
    mult = format["multiply"](["", ""])
    grouping = format["grouping"]
    group = grouping("")
    format_G = format["geometry tensor"]
    gl = group[0]
    gr = group[1]

    # Get the number of geometry declaration, possibly offset value
    num_geo = offset + len(geo_terms)
    new_prods = []

    # Split the expression into products
    prods = split_expression(expression, format, add)
    consts = []

    # Loop products and check if the variables are constant
    for p in prods:
        vrs = split_expression(p, format, mult)
        geos = []

        # Generate geo code for constant coefficients e.g., w[0][5]
        new_vrs = []
        for v in vrs:

            # If variable is a group, get the geometry terms and
            # update geo number
            if v[0] == gl and v[-1] == gr:
                v = get_geo_terms(v[1:-1], geo_terms, offset, format)
                num_geo = offset + len(geo_terms)

                # If we still have a sum, regroup
                if len(v.split(add)) > 1:
                    v = grouping(v)

            # Append to new variables
            new_vrs.append(v)

            # If variable is constants, add to geo terms
            constant = is_constant(v, format)
            if constant:
                geos.append(v)

        # Update variable list
        vrs = new_vrs
        vrs.sort()

        # Sort geo and create geometry term
        geos.sort()
        geo = mult.join(geos)

        # Handle geometry term appropriately
        if geo:
            if geos != vrs:
                if len(geos) > 1:
                    for g in geos:
                        vrs.remove(g)
                    if geo not in geo_terms:
                        geo_terms[geo] = format_G + str(num_geo)
                        num_geo += 1
                    vrs.append(geo_terms[geo])
                new_prods.append(mult.join(vrs))
            else:
                consts.append(mult.join(vrs))
        else:
            new_prods.append(mult.join(vrs))

    if consts:
        if len(consts) > 1:
            c = grouping(add.join(consts))
        else:
            c = add.join(consts)
        if c not in geo_terms:
            geo_terms[c] = format_G + str(num_geo)
            num_geo += 1
        consts = [geo_terms[c]]

    return add.join(new_prods + consts)


def get_constants(expression, const_terms, format, constants=[]):
    """This function returns a new expression where all geometry terms
    have been substituted with geometry declarations, these
    declarations are added to the const_terms dictionary.

    """

    # Get formats
    add = format["add"](["", ""])
    mult = format["multiply"](["", ""])
    grouping = format["grouping"]
    format_G = format["geometry tensor"] + "".join(constants)  # format["geometry tensor"]

    # Get the number of geometry declaration, possibly offset value
    num_geo = len(const_terms)
    new_prods = []

    # Split the expression into products
    prods = split_expression(expression, format, add)
    consts = []

    # Loop products and check if the variables are constant
    for p in prods:
        vrs = split_expression(p, format, mult)
        geos = []

        # Generate geo code for constant coefficients e.g., w[0][5]
        new_vrs = []
        for v in vrs:

            # If variable is constants, add to geo terms
            constant = is_constant(v, format, constants)
            if constant:
                geos.append(v)
            # Append to new variables
            new_vrs.append(v)

        # Update variable list
        vrs = new_vrs
        vrs.sort()

        # Sort geo and create geometry term
        geos.sort()
        geo = mult.join(geos)
        if geo:
            if geos != vrs:
                for g in geos:
                    vrs.remove(g)
                if geo not in const_terms:
                    const_terms[geo] = format_G + str(num_geo)
                    num_geo += 1
                vrs.append(const_terms[geo])
                new_prods.append(mult.join(vrs))
            else:
                consts.append(mult.join(vrs))
        else:
            new_prods.append(mult.join(vrs))

    if consts:
        if len(consts) > 1:
            c = grouping(add.join(consts))
        else:
            c = add.join(consts)
        if c not in const_terms:
            const_terms[c] = format_G + str(num_geo)
            num_geo += 1
        consts = [const_terms[c]]

    return add.join(new_prods + consts)


def get_indices(variable, format, from_get_indices=False):
    """This function returns the indices of a given variable. E.g.,
    P[0][j],            returns ['j']
    P[ip][k],           returns ['ip','k']
    P[ip][nzc0[j] + 3], returns ['ip','j']
    w[0][j + 2]         , returns [j]"""

    add = format["add"](["", ""])
    mult = format["multiply"](["", ""])
    format_access = format["array access"]
    access = format_access("")

    l = access[0]  # noqa: E741
    r = access[1]

    indices = []

    # If there are no '[' in variable and self is the caller
    if not variable.count(l) and from_get_indices:
        adds = split_expression(variable, format, add)
        for a in adds:
            mults = split_expression(a, format, mult)
            for m in mults:
                try:
                    float(m)
                except Exception:
                    if m not in indices:
                        indices.append(m)
    else:
        index = ""
        left = 0
        inside = False
        # Loop all characters in variable and find indices
        for c in variable:
            # If character is ']' reduce left count
            if c == r:
                left -= 1

            # If the '[' count has returned to zero, we have a
            # complete index
            if left == 0 and inside:
                try:
                    eval(index)
                except Exception:
                    indices += get_indices(index, format, True)
                index = ""
                inside = False

            # If we're inside an access, add character to index
            if inside:
                index += c

            # If character is '[' increase the count, and we're inside
            # an access
            if c == l:
                inside = True
                left += 1

    return indices


def get_variables(expression, variables, format, constants=[]):
    """This function returns a new expression where all geometry terms
    have been substituted with geometry declarations, these
    declarations are added to the const_terms dictionary.

    """

    # Get formats
    add = format["add"](["", ""])
    mult = format["multiply"](["", ""])
    format_access = format["array access"]
    access = format_access("")
    format_F = format["function value"]

    l = access[0]  # noqa: E741

    # If we don't have any access operators in expression,
    # we don't have any variables
    if expression.count(l) == 0:
        return expression

    # Get the number of geometry declaration, possibly offset value
    num_var = len(variables)
    new_prods = []
    used_vars = []

    # Split the expression into products
    prods = split_expression(expression, format, add)

    # Loop products and check if the variables are constant
    for p in prods:
        vrs = split_expression(p, format, mult)
        # Variables with respect to the constants in list
        variables_of_interest = []

        # Generate geo code for constant coefficients e.g., w[0][5]
        new_vrs = []
        for v in vrs:
            # If we don't have any access operators, we don't have a
            # variable
            if v.count(l) == 0:
                new_vrs.append(v)
                continue

            # Check if we have a variable that depends on one of the
            # constants First check the easy way
            is_var = False
            for c in constants:
                if format_access(c) in v:
                    is_var = True
                    break
            if is_var:
                variables_of_interest.append(v)
                continue

            # Then check the hard way
            # Get list of indices
            indices = get_indices(v, format)
            depends = [True for c in constants if c in indices]
            if any(depends):
                variables_of_interest.append(v)
            else:
                new_vrs.append(v)

        variables_of_interest.sort()
        variables_of_interest = mult.join(variables_of_interest)

        # If we have some variables, declare new variable if needed
        # and add to list of variables
        if variables_of_interest:
            # If we didn't already declare this variable do so
            if variables_of_interest not in variables:
                variables[variables_of_interest] = format_F + str(num_var)
                num_var += 1

            # Get mapped variable
            mv = variables[variables_of_interest]
            new_vrs.append(mv)
            if mv not in used_vars:
                used_vars.append(mv)

        # Sort variables and add to list of products
        new_vrs.sort()
        new_prods.append(mult.join(new_vrs))

    # Sort list of products and return the sum
    new_prods.sort()
    return (add.join(new_prods), used_vars)
