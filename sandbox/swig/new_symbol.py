"Some simple functions for manipulating expressions symbolically"

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2009-07-12 -- 2009-07-15"
__copyright__ = "Copyright (C) 2009 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# FFC common modules
from ffc.common.log import debug, error

#import psyco
#psyco.full()

# TODO: Use proper errors, not just RuntimeError
# TODO: Change all if value == 0.0 to something more safe

# Some basic variables
BASIS = 0
IP  = 1
GEO = 2
CONST = 3
type_to_string = {BASIS:"BASIS", IP:"IP",GEO:"GEO", CONST:"CONST"}
format = None

# Functions and dictionaries for cache implementation
# Increases speed and should also reduce memory consumption
_float_cache = {}
def create_float(val):
    if val in _float_cache:
#        print "found %f in cache" %val
        return _float_cache[val]
    float_val = FloatValue(val)
    _float_cache[val] = float_val
    return float_val

_symbol_cache = {}
def create_symbol(variable, symbol_type, base_expr=None, base_op=0):
    key = (variable, symbol_type, base_expr, base_op)
    if key in _symbol_cache:
#        print "found %s in cache" %variable
        return _symbol_cache[key]
    symbol = Symbol(variable, symbol_type, base_expr, base_op)
    _symbol_cache[key] = symbol
    return symbol

_product_cache = {}
def create_product(variables):
#    variables.sort()
#    key = tuple(sorted(variables))
    key = tuple(variables)
    if key in _product_cache:
#        print "found %s in cache" %str(key)
        return _product_cache[key]
    product = Product(key)
    _product_cache[key] = product
    return product

_sum_cache = {}
def create_sum(variables):
#    variables.sort()
#    key = tuple(sorted(variables))
    key = tuple(variables)
    if key in _sum_cache:
#        print "found %s in cache" %str(key)
        return _sum_cache[key]
    s = Sum(key)
    _sum_cache[key] = s
    return s

_fraction_cache = {}
def create_fraction(num, denom):
    key = (num, denom)
    if key in _fraction_cache:
#        print "found %s in cache" %str(key)
        return _fraction_cache[key]
    fraction = Fraction(num, denom)
    _fraction_cache[key] = fraction
    return fraction

# Function to set global format to avoid passing around the dictionary
def set_format(_format):
    global format
    format = _format
    set_format_float(format)
    set_format_prod(format)
    set_format_sum(format)
    set_format_frac(format)
    global format_comment
    format_comment = format["comment"]

def get_format():
    return format

def generate_aux_constants(constant_decl, name, var_type, print_ops=False):
    "A helper tool to generate code for constant declarations"
    code = []
    append = code.append
    ops = 0
    sorted_list = sorted([(v, k) for k, v in constant_decl.iteritems()])
    for s in sorted_list:
        c = s[1]
#        print "k: ", s[0]
#        print "c: ", repr(c)
#        debug("c orig: " + str(c))
        c = c.expand().reduce_ops()
#        print "ops: ", c.ops()
#        debug("c opt:  " + str(c))
        if print_ops:
            op = c.ops()
            ops += op
#            append(format["comment"]("Number of operations: %d" %c.ops()))
            append(format_comment("Number of operations: %d" %op))
            append((var_type + name + str(s[0]), str(c)))
            append("")
        else:
            ops += c.ops()
            append((var_type + name + str(s[0]), str(c)))

    return (ops, code)

def optimise_code(expr, ip_consts, geo_consts, trans_set):
    """Optimise a given expression with respect to, basis functions,
    integration points variables and geometric constants.
    The function will update the dictionaries ip_const and geo_consts with new
    declarations and update the trans_set (used transformations)."""

    format_G  = format["geometry tensor"]
    format_ip = format["integration points"]

    trans_set_update = trans_set.update

    # Return constant symbol if value is zero
    if expr.val == 0.0:
#        return FloatValue(0)
        return create_float(0)

    # Reduce expression with respect to basis function variable
#    debug("\n\nexpr before exp: " + repr(expr))
#    expr = expr.expand()
#    debug("\n\nexpr: " + str(expr))
#    debug("\n\nexpr: " + repr(expr))
#    basis_expressions = expr.reduce_vartype(BASIS)
    basis_expressions = expr.expand().reduce_vartype(BASIS)

    # If we had a product instance we'll get a tuple back so embed in list
    if not isinstance(basis_expressions, list):
        basis_expressions = [basis_expressions]

    basis_vals = []
#    basis_vals_append = basis_vals.append
    # Process each instance of basis functions
#    for b in basis_expressions:
    for basis, ip_expr in basis_expressions:
        # Get the basis and the ip expression
#        basis, ip_expr = b
#        debug("\nbasis\n" + str(basis))
#        debug("ip_epxr\n" + str(ip_expr))

        # If we have no basis (like functionals) create a const
        if not basis:
#            basis = FloatValue(1)
            basis = create_float(1)

#        if Product([basis, ip_expr], False).expand() != expr:

#            prod = Product([basis, ip_expr], False).expand()
#            print "prod == sum: ", isinstance(prod, Sum)
#            print "expr == sum: ", isinstance(expr, Sum)

#            print "prod.pos: ", prod.pos
#            print "expr.pos: ", expr.pos
#            print "expr.pos = prod.pos: ", expr.pos == prod.pos

#            print "prod.neg: ", prod.neg
#            print "expr.neg: ", expr.neg
#            print "expr.neg = prod.neg: ", expr.neg == prod.neg

#            print "equal: ", prod == expr

#            print "\nprod:    ", prod
#            print "\nexpr:    ", expr
#            print "\nbasis:   ", basis
#            print "\nip_expr: ", ip_expr
#            raise RuntimeError("Not equal")

        # If the ip expression doesn't contain any operations skip remainder
        if not ip_expr:
            basis_vals.append(basis)
#            basis_vals_append(basis)
            continue
        if not ip_expr.ops() > 0:
#            basis_vals.append(Product([basis, ip_expr]))
#            basis_vals_append(Product([basis, ip_expr]))
            basis_vals.append(create_product([basis, ip_expr]))
#            basis_vals_append(create_product([basis, ip_expr]))
            continue

        # Reduce the ip expressions with respect to IP variables
#        ip_expr = ip_expr.expand()
#        ip_expressions = ip_expr.reduce_vartype(IP)
        ip_expressions = ip_expr.expand().reduce_vartype(IP)

        # If we had a product instance we'll get a tuple back so embed in list
        if not isinstance(ip_expressions, list):
            ip_expressions = [ip_expressions]

        ip_vals = []
#        ip_vals_append = ip_vals.append
        # Loop ip expressions
        for ip in ip_expressions:
            ip_dec, geo = ip
#            debug("\nip_dec: " + str(ip_dec))
#            debug("\ngeo: " + str(geo))
            # Update transformation set with those values that might be
            # embedded in IP terms
            if ip_dec:
#                trans_set.update(map(lambda x: str(x), ip_dec.get_unique_vars(GEO)))
                trans_set_update(map(lambda x: str(x), ip_dec.get_unique_vars(GEO)))

            # Append and continue if we did not have any geo values
            if not geo:
                ip_vals.append(ip_dec)
#                ip_vals_append(ip_dec)
                continue

            # Update the transformation set with the variables in the geo term
#            trans_set.update(map(lambda x: str(x), geo.get_unique_vars(GEO)))
            trans_set_update(map(lambda x: str(x), geo.get_unique_vars(GEO)))

            # Only declare auxiliary geo terms if we can save operations            
            if geo.ops() > 0:
#                debug("geo: " + str(geo))
                # If the geo term is not in the dictionary append it
                if not geo_consts.has_key(geo):
                    geo_consts[geo] = len(geo_consts)

                # Substitute geometry expression
#                geo = Symbol(format_G + str(geo_consts[geo]), GEO)
                geo = create_symbol(format_G + str(geo_consts[geo]), GEO)

            # If we did not have any ip_declarations use geo, else create a
            # product and append to the list of ip_values
            if not ip_dec:
                ip_dec = geo
            else:
#                ip_dec = Product([ip_dec, geo])
                ip_dec = create_product([ip_dec, geo])
            ip_vals.append(ip_dec)
#            ip_vals_append(ip_dec)

        # Create sum of ip expressions to multiply by basis
        if len(ip_vals) > 1:
#            ip_expr = Sum(ip_vals)
            ip_expr = create_sum(ip_vals)
        elif ip_vals:
            ip_expr = ip_vals.pop()

        # If we can save operations by declaring it as a constant do so, if it
        # is not in IP dictionary, add it and use new name
        if ip_expr.ops() > 0:
            if not ip_expr in ip_consts:
                ip_consts[ip_expr] = len(ip_consts)

            # Substitute ip expression
#            ip_expr = Symbol(format_G + format_ip + str(ip_consts[ip_expr]), IP)
            ip_expr = create_symbol(format_G + format_ip + str(ip_consts[ip_expr]), IP)

        # Multiply by basis and append to basis vals
#        basis_vals.append(Product([basis, ip_expr]).expand())
#        basis_vals_append(Product([basis, ip_expr]).expand())
        basis_vals.append(create_product([basis, ip_expr]).expand())
#        basis_vals_append(create_product([basis, ip_expr]).expand())

    # Return sum of basis values
#    return Sum(basis_vals)
    return create_sum(basis_vals)

#def reduce_vartype(o, var_type):

#    if o._class in ("float", "sym"):
#        if o.t == var_type:
##            return (self, FloatValue(1))
#            return (o, create_float(1))
#        return ((), o)

#    elif o._class  == "prod":
#        # Sort variables according to type
#        found = []
#        found_append = found.append
#        remains = []
#        remains_append = remains.append
#        for v in o.vrs:
#            if v.t == var_type:
#                found_append(v)
#                continue
#            remains_append(v)

#        # Create appropriate object for found
#        if len(found) > 1:
#            found = create_product(found)
#        elif found:
#            found = found.pop()
#        # We did not find any variables
#        else:
#            return ((), o)
#        # Create appropriate object for remains
#        if len(remains) > 1:
#            remains = create_product(remains)
#        elif remains:
#            remains = remains.pop()
#        # We don't have anything left
#        else:
#            return (o, create_float(1))
#        # Return whatever we found
#        return (found, remains)

#    elif o._class  == "sum":
#        found = {}
#        # Loop members and reduce them by vartype
#        for v in o.pos + o.neg:
#            f, r = reduce_vartype(v, var_type)
#            if f in found:
#                found[f].append(r)
#            else:
#                found[f] = [r]

#        # Create the return value
#        returns = []
#        append = returns.append
#        for f, r in found.iteritems():
#            if len(r) > 1:
#                # Use expand to group expressions
#                r = create_sum(r).expand()
#            elif r:
#                r = r.pop()
#            append((f, r))
#        return returns

#    # Fraction
#    num_found, num_remain = reduce_vartype(o.num, var_type)

#    # TODO: Remove this test later, expansion should have taken care of
#    # no denominator
#    if not o.denom:
#        raise RuntimeError("This fraction should have been expanded")

#    # If the denominator is not a Sum things are straightforward
#    denom_found = None
#    denom_remain = None
##        if not isinstance(self.denom, Sum):
#    if o.denom._class != "sum":
#        denom_found, denom_remain = reduce_vartype(o.denom, var_type)

#    # If we have a Sum in the denominator, all terms must be reduced by
#    # the same terms to make sense
#    else:
#        remain = []
#        for m in o.denom.pos + o.denom.neg:
#            d_found, d_remain = reduce_vartype(m, var_type)
#            # If we've found a denom, but the new found is different from
#            # the one already found, terminate loop since it wouldn't make
#            # sense to reduce the fraction
#            if denom_found != None and str(d_found) != str(denom_found):
#                # In case we did not find any variables of given type in the numerator
#                # declare a constant. We always have a remainder.
#                return (num_found, create_fraction(num_remain, o.denom))

#            denom_found = d_found
#            remain.append(d_remain)

#        # There is always a non-const remainder if denominator was a sum
#        denom_remain = create_sum(remain)

#    # If we have found a common denominator, but no found numerator,
#    # create a constant
#    # TODO: Add more checks to avoid expansion
#    found = None
#    # There is always a remainder
#    remain = create_fraction(num_remain, denom_remain).expand()

#    if num_found:
#        if denom_found:
#            found = create_fraction(num_found, denom_found)
#        else:
#            found = num_found
#    else:
#        if denom_found:
#            found = create_fraction(create_float(1), denom_found)
#        else:
#            found = ()
#    return (found, remain)

from floatvalue_obj import FloatValue, set_format as set_format_float
from symbol_obj     import Symbol
from product_obj    import Product, set_format as set_format_prod
from sum_obj        import Sum, group_fractions, set_format as set_format_sum
from fraction_obj   import Fraction, set_format as set_format_frac


