# -*- coding: utf-8 -*-
"This file contains functions to optimise the code generated for quadrature representation."

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

# FFC modules
from ffc.log import error
from ffc.quadrature.cpp import format

# TODO: Use proper errors, not just RuntimeError.
# TODO: Change all if value == 0.0 to something more safe.

# Some basic variables.
BASIS = 0
IP = 1
GEO = 2
CONST = 3
type_to_string = {BASIS: "BASIS", IP: "IP", GEO: "GEO", CONST: "CONST"}

# Functions and dictionaries for cache implementation.
# Increases speed and should also reduce memory consumption.
_float_cache = {}


def create_float(val):
    if val in _float_cache:
        return _float_cache[val]
    float_val = FloatValue(val)
    _float_cache[val] = float_val
    return float_val


_symbol_cache = {}


def create_symbol(variable, symbol_type, base_expr=None, base_op=0):
    key = (variable, symbol_type, base_expr, base_op)
    if key in _symbol_cache:
        return _symbol_cache[key]
    symbol = Symbol(variable, symbol_type, base_expr, base_op)
    _symbol_cache[key] = symbol
    return symbol


_product_cache = {}


def create_product(variables):
    # NOTE: If I switch on the sorted line, it might be possible to find more
    # variables in the cache, but it adds some overhead so I don't think it
    # pays off. The member variables are also sorted in the classes
    # (Product and Sum) so the list 'variables' is probably already sorted.
    key = tuple(variables)
    if key in _product_cache:
        return _product_cache[key]
    product = Product(key)
    _product_cache[key] = product
    return product


_sum_cache = {}


def create_sum(variables):
    # NOTE: If I switch on the sorted line, it might be possible to
    # find more variables in the cache, but it adds some overhead so I
    # don't think it pays off. The member variables are also sorted in
    # the classes (Product and Sum) so the list 'variables' is
    # probably already sorted.
    key = tuple(variables)
    if key in _sum_cache:
        return _sum_cache[key]
    s = Sum(key)
    _sum_cache[key] = s
    return s


_fraction_cache = {}


def create_fraction(num, denom):
    key = (num, denom)
    if key in _fraction_cache:
        return _fraction_cache[key]
    fraction = Fraction(num, denom)
    _fraction_cache[key] = fraction
    return fraction


def generate_aux_constants(constant_decl, name, var_type, print_ops=False):
    "A helper tool to generate code for constant declarations."
    format_comment = format["comment"]
    code = []
    append = code.append
    ops = 0
    for num, expr in sorted((v, k) for k, v in sorted_by_key(constant_decl)):
        # Expand and reduce expression (If we don't already get
        # reduced expressions.)
        expr = expr.expand().reduce_ops()
        if print_ops:
            op = expr.ops()
            ops += op
            append(format_comment("Number of operations: %d" % op))
            append(var_type(name(num), str(expr)))
            append("")
        else:
            ops += expr.ops()
            append(var_type(name(num), str(expr)))

    return (ops, code)


def optimise_code(expr, ip_consts, geo_consts, trans_set):
    """Optimise a given expression with respect to, basis functions,
    integration points variables and geometric constants.  The
    function will update the dictionaries ip_const and geo_consts with
    new declarations and update the trans_set (used
    transformations).

    """

    format_G = format["geometry constant"]
    format_I = format["ip constant"]
    trans_set_update = trans_set.update

    # Return constant symbol if expanded value is zero.
    exp_expr = expr.expand()
    if exp_expr.val == 0.0:
        return create_float(0)

    # Reduce expression with respect to basis function variable.
    basis_expressions = exp_expr.reduce_vartype(BASIS)

    # If we had a product instance we'll get a tuple back so embed in
    # list.
    if not isinstance(basis_expressions, list):
        basis_expressions = [basis_expressions]

    basis_vals = []
    # Process each instance of basis functions.
    for basis, ip_expr in basis_expressions:
        # Get the basis and the ip expression.

        # If we have no basis (like functionals) create a const.
        if not basis:
            basis = create_float(1)

        # If the ip expression doesn't contain any operations skip
        # remainder
        if not ip_expr or ip_expr.val == 0.0:
            basis_vals.append(basis)
            continue
        if not ip_expr.ops() > 0:
            basis_vals.append(create_product([basis, ip_expr]))
            continue

        # Reduce the ip expressions with respect to IP variables.
        ip_expressions = ip_expr.expand().reduce_vartype(IP)

        # If we had a product instance we'll get a tuple back so embed in list.
        if not isinstance(ip_expressions, list):
            ip_expressions = [ip_expressions]

        ip_vals = []
        # Loop ip expressions.
        for ip in sorted(ip_expressions):
            ip_dec, geo = ip

            # Update transformation set with those values that might
            # be embedded in IP terms.
            if ip_dec and ip_dec.val != 0.0:
                trans_set_update([str(x) for x in ip_dec.get_unique_vars(GEO)])

            # Append and continue if we did not have any geo values.
            if not geo or geo.val == 0.0:
                if ip_dec and ip_dec.val != 0.0:
                    ip_vals.append(ip_dec)
                continue

            # Update the transformation set with the variables in the
            # geo term.
            trans_set_update([str(x) for x in geo.get_unique_vars(GEO)])

            # Only declare auxiliary geo terms if we can save
            # operations.
            if geo.ops() > 0:
                # If the geo term is not in the dictionary append it.
                if geo not in geo_consts:
                    geo_consts[geo] = len(geo_consts)

                # Substitute geometry expression.
                geo = create_symbol(format_G(geo_consts[geo]), GEO)

            # If we did not have any ip_declarations use geo, else
            # create a product and append to the list of ip_values.
            if not ip_dec or ip_dec.val == 0.0:
                ip_dec = geo
            else:
                ip_dec = create_product([ip_dec, geo])
            ip_vals.append(ip_dec)

        # Create sum of ip expressions to multiply by basis.
        if len(ip_vals) > 1:
            ip_expr = create_sum(ip_vals)
        elif ip_vals:
            ip_expr = ip_vals.pop()

        # If we can save operations by declaring it as a constant do
        # so, if it is not in IP dictionary, add it and use new name.
        if ip_expr.ops() > 0 and ip_expr.val != 0.0:
            if ip_expr not in ip_consts:
                ip_consts[ip_expr] = len(ip_consts)

            # Substitute ip expression.
            ip_expr = create_symbol(format_I(ip_consts[ip_expr]), IP)

        # Multiply by basis and append to basis vals.
        basis_vals.append(create_product([basis, ip_expr]))

    # Return (possible) sum of basis values.
    if len(basis_vals) > 1:
        return create_sum(basis_vals)
    elif basis_vals:
        return basis_vals[0]

    # Where did the values go?
    error("Values disappeared.")


from .floatvalue import FloatValue
from .symbol import Symbol
from .product import Product
from .sumobj import Sum
from .fraction import Fraction
