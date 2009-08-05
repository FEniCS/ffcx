"Some simple functions for manipulating expressions symbolically"

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2009-07-12 -- 2009-07-15"
__copyright__ = "Copyright (C) 2009 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# FFC common modules
from ffc.common.log import debug, error

BASIS = 0
IP  = 1
GEO = 2
CONST = 3
type_to_string = {BASIS:"BASIS", IP:"IP",GEO:"GEO", CONST:"CONST"}

format = None
EPS = 1e-12

# TODO: Use proper errors, not just RuntimeError
# TODO: Change all if value == 0.0 to something more safe
from floatvalue_obj import FloatValue, set_format as set_format_float
from symbol_obj     import Symbol
from product_obj    import Product, set_format as set_format_prod
from sum_obj        import Sum, group_fractions, set_format as set_format_sum
from fraction_obj   import Fraction, set_format as set_format_frac

def set_format(_format):
    global format
    format = _format
    set_format_float(format)
    set_format_prod(format)
    set_format_sum(format)
    set_format_frac(format)

def get_format():
    return format


def generate_aux_constants(constant_decl, name, var_type, print_ops=False):
    "A helper tool to generate code for constant declarations"
    code = []
    append = code.append
    ops = 0
    sorted_list = [(v, k) for k, v in constant_decl.iteritems()]
    sorted_list.sort()
    for s in sorted_list:
        c = s[1]
#        debug("c orig: " + str(c))
#        c = c.expand().reduce_ops()
#        debug("c opt:  " + str(c))
        ops += c.ops()
        if print_ops:
            append(format["comment"]("Number of operations: %d" %c.ops()))
            append((var_type + name + str(s[0]), str(c)))
            append("")
        else:
            append((var_type + name + str(s[0]), str(c)))
    return (ops, code)

def optimise_code(expr, ip_consts, geo_consts, trans_set):
    """Optimise a given expression with respect to, basis functions,
    integration points variables and geometric constants.
    The function will update the dictionaries ip_const and geo_consts with new
    declarations and update the trans_set (used transformations)."""

    format_G  = format["geometry tensor"]
    format_ip = format["integration points"]

    # Return constant symbol if value is zero
    if expr.val == 0.0:
        return FloatValue(0)

    # Reduce expression with respect to basis function variable
    debug("\n\nexpr before exp: " + repr(expr))
    expr = expr.expand()
    debug("\n\nexpr: " + str(expr))
    debug("\n\nexpr: " + repr(expr))
    basis_expressions = expr.reduce_vartype(BASIS)

    # If we had a product instance we'll get a tuple back so embed in list
    if not isinstance(basis_expressions, list):
        basis_expressions = [basis_expressions]

    basis_vals = []
    # Process each instance of basis functions
    for b in basis_expressions:
        # Get the basis and the ip expression
        basis, ip_expr = b
        debug("\nbasis\n" + str(basis))
        debug("ip_epxr\n" + str(ip_expr))

        # If we have no basis (like functionals) create a const
        if not basis:
            basis = FloatValue(1)

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
            continue
        if not ip_expr.ops() > 0:
            basis_vals.append(Product([basis, ip_expr]))
            continue

        # Reduce the ip expressions with respect to IP variables
        ip_expr = ip_expr.expand()
        ip_expressions = ip_expr.reduce_vartype(IP)

        # If we had a product instance we'll get a tuple back so embed in list
        if not isinstance(ip_expressions, list):
            ip_expressions = [ip_expressions]

        ip_vals = []
        # Loop ip expressions
        for ip in ip_expressions:
            ip_dec, geo = ip
            debug("\nip_dec: " + str(ip_dec))
            debug("\ngeo: " + str(geo))
            # Update transformation set with those values that might be
            # embedded in IP terms
            if ip_dec:
                trans_set.update(map(lambda x: str(x), ip_dec.get_unique_vars(GEO)))

            # Append and continue if we did not have any geo values
            if not geo:
                ip_vals.append(ip_dec)
                continue

            # Update the transformation set with the variables in the geo term
            trans_set.update(map(lambda x: str(x), geo.get_unique_vars(GEO)))

            # Only declare auxiliary geo terms if we can save operations            
            if geo.ops() > 0:
                debug("geo: " + str(geo))
                # If the geo term is not in the dictionary append it
                if not geo_consts.has_key(geo):
                    geo_consts[geo] = len(geo_consts)

                # Substitute geometry expression
                geo = Symbol(format_G + str(geo_consts[geo]), GEO)

            # If we did not have any ip_declarations use geo, else create a
            # product and append to the list of ip_values
            if not ip_dec:
                ip_dec = geo
            else:
                ip_dec = Product([ip_dec, geo])
            ip_vals.append(ip_dec)

        # Create sum of ip expressions to multiply by basis
        if len(ip_vals) > 1:
            ip_expr = Sum(ip_vals)
        elif ip_vals:
            ip_expr = ip_vals.pop()

        # If we can save operations by declaring it as a constant do so, if it
        # is not in IP dictionary, add it and use new name
        if ip_expr.ops() > 0:
            if not ip_expr in ip_consts:
                ip_consts[ip_expr] = len(ip_consts)

            # Substitute ip expression
            ip_expr = Symbol(format_G + format_ip + str(ip_consts[ip_expr]), IP)

        # Multiply by basis and append to basis vals
        basis_vals.append(Product([basis, ip_expr]).expand())

    # Return sum of basis values
    return Sum(basis_vals)



