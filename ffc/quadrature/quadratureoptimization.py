__author__ = "Kristian B. Oelgaard (k.b.oelgaard@gmail.com)"
__date__ = "2010-02-08"
__copyright__ = "Copyright (C) 2010 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# Last changed: 2010-03-11

# FFC modules
from ffc.log import info
from ffc.cpp import format
from ffc.quadrature.symbolics import optimise_code, BASIS, IP, GEO, CONST
from ffc.quadrature.symbolics import create_product, create_sum, create_symbol, create_fraction

def optimize_integral_ir(ir):
    "Compute optimized intermediate representation of integral."

    parameters = ir["optimise_parameters"]
    if parameters["optimisation"]:
        integrals  = ir["trans_integrals"]
        domain_type =ir["domain_type"]
        num_facets = ir["num_facets"]
        geo_consts = ir["geo_consts"]
        if domain_type == "cell":
            info("Optimising expressions for cell integral")
            if parameters["optimisation"] in ("precompute_ip_const", "precompute_basis_const"):
                _precompute_expressions(integrals, geo_consts, parameters["optimisation"])
            else:
                _simplify_expression(integrals, geo_consts)
        elif domain_type == "exterior_facet":
            for i in range(num_facets):
                info("Optimising expressions for facet integral %d" % i)
                if parameters["optimisation"] in ("precompute_ip_const", "precompute_basis_const"):
                    _precompute_expressions(integrals[i], geo_consts, parameters["optimisation"])
                else:
                    _simplify_expression(integrals[i], geo_consts)
        elif domain_type == "interior_facet":
            for i in range(num_facets):
                for j in range(num_facets):
                    info("Optimising expressions for facet integral (%d, %d)" % (i, j))
                    if parameters["optimisation"] in ("precompute_ip_const", "precompute_basis_const"):
                        _precompute_expressions(integrals[i][j], geo_consts,parameters["optimisation"])
                    else:
                        _simplify_expression(integrals[i][j], geo_consts)
        else:
            error("Unhandled domain type: " + str(domain_type))

    return ir

def _simplify_expression(integral, geo_consts):
    for points, terms, functions, ip_consts, coordinate in integral:
        for loop, (data, entry_vals) in terms.iteritems():
            t_set, u_weights, u_psi_tables, u_nzcs, basis_consts = data
            new_entry_vals = []
            for entry, val, ops in entry_vals:
                value = optimise_code(val, ip_consts, geo_consts, t_set)
                # Check if value is zero
                if value.val:
                    new_entry_vals.append((entry, value, value.ops()))
            terms[loop][1] = new_entry_vals

def _precompute_expressions(integral, geo_consts, optimisation):
    for points, terms, functions, ip_consts, coordinate in integral:
        for loop, (data, entry_vals) in terms.iteritems():
            t_set, u_weights, u_psi_tables, u_nzcs, basis_consts = data
            new_entry_vals = []
            for entry, val, ops in entry_vals:
                value = _extract_variables(val, basis_consts, ip_consts, geo_consts, t_set, optimisation)
                # Check if value is zero
                if value.val:
                    new_entry_vals.append((entry, value, value.ops()))
            terms[loop][1] = new_entry_vals

def _extract_variables(val, basis_consts, ip_consts, geo_consts, t_set, optimisation):
    f_G = format["geometry constant"]
    f_I = format["ip constant"]
    f_B = format["basis constant"]

    if val._prec == 0:
        return val
    elif val._prec == 1:
        if val.base_expr is None:
            return val
        new_base = _extract_variables(val.base_expr, basis_consts, ip_consts, geo_consts, t_set, optimisation)
        new_sym = create_symbol(val.v, val.t, new_base, val.base_op)
        new_sym.exp = val.exp
        if new_sym.t == BASIS:
            return _reduce_expression(new_sym, [], basis_consts, f_B, True)
        elif new_sym.t == IP:
            return _reduce_expression(new_sym, [], ip_consts, f_I, True)
        elif new_sym.t == GEO:
            return _reduce_expression(new_sym, [], geo_consts, f_G, True)
    # First handle child classes of product and sum.
    elif val._prec in (2, 3):
        new_vars = []
        for v in val.vrs:
            new_vars.append(_extract_variables(v, basis_consts, ip_consts, geo_consts, t_set, optimisation))
        if val._prec == 2:
            new_val = create_product(new_vars)
        if val._prec == 3:
            new_val = create_sum(new_vars)
    elif val._prec == 4:
        num = _extract_variables(val.num, basis_consts, ip_consts, geo_consts, t_set, optimisation)
        denom = _extract_variables(val.denom, basis_consts, ip_consts, geo_consts, t_set, optimisation)
        return create_fraction(num, denom)
    else:
        error("Unknown symbolic type: %s" % repr(val))

    # Sort variables of product and sum.
    b_c, i_c, g_c = [], [], []
    for v in new_val.vrs:
        if v.t == BASIS:
            if optimisation == "precompute_basis_const":
                b_c.append(v)
        elif v.t == IP:
            i_c.append(v)
        else:
            g_c.append(v)
    vrs = new_val.vrs[:]
    for v in g_c + i_c + b_c:
        vrs.remove(v)
    i_c.extend(_reduce_expression(new_val, g_c, geo_consts, f_G))
    vrs.extend(_reduce_expression(new_val, i_c, ip_consts, f_I))
    vrs.extend(_reduce_expression(new_val, b_c, basis_consts, f_B))

#    print "b_c: "
#    for b in b_c:
#        print b
#    print "basis"
#    for k,v in basis_consts.items():
#        print "k: ", k
#        print "v: ", v
#    print "geo"
#    for k,v in geo_consts.items():
#        print "k: ", k
#        print "v: ", v
#    print "ret val: ", val

    if len(vrs) > 1:
        if new_val._prec == 2:
            new_object = create_product(vrs)
        elif new_val._prec == 3:
            new_object = create_sum(vrs)
        else:
            error("Must have product or sum here: %s" % repr(new_val))
        if new_object.t == BASIS:
            if optimisation == "precompute_ip_const":
                return new_object
            elif optimisation == "precompute_basis_const":
                return _reduce_expression(new_object, [], basis_consts, f_B, True)
        elif new_object.t == IP:
            return _reduce_expression(new_object, [], ip_consts, f_I, True)
        elif new_object.t == GEO:
            return _reduce_expression(new_object, [], geo_consts, f_G, True)
    return vrs[0]

#    if new_val._prec == 2:
#        if len(vrs) > 1:
#            new_prod = create_product(vrs)
#            if new_prod.t == BASIS:
#                if optimisation == "precompute_ip_const":
#                    return new_prod
#                elif optimisation == "precompute_basis_const":
#                    return _reduce_expression(new_prod, [], basis_consts, f_B, True)
#            elif new_prod.t == IP:
#                return _reduce_expression(new_prod, [], ip_consts, f_I, True)
#            elif new_prod.t == GEO:
#                return _reduce_expression(new_prod, [], geo_consts, f_G, True)
#        return vrs[0]
#    elif new_val._prec == 3:
#        if len(vrs) > 1:
#            new_sum = create_sum(vrs)
#            if new_sum.t == BASIS:
#                return new_sum
##                return _reduce_expression(new_sum, [], basis_consts, f_B, True)
#            elif new_sum.t == IP:
#                return _reduce_expression(new_sum, [], ip_consts, f_I, True)
#            elif new_sum.t == GEO:
#                return _reduce_expression(new_sum, [], geo_consts, f_G, True)
#        return vrs[0]
#    else:
#        error("Must have product or sum here: %s" % repr(new_val))

def _reduce_expression(expr, symbols, const_dict, f_name, use_expr_type=False):
    if use_expr_type:
        if expr not in const_dict:
            const_dict[expr] = len(const_dict)
        return create_symbol(f_name(const_dict[expr]), expr.t)
    # Only something to be done if we have more than one symbol.
    if len(symbols) > 1:
        sym_type = symbols[0].t
        # Create new symbol.
        if expr._prec == 2:
            new_sym = create_product(symbols)
        elif expr._prec == 3:
            new_sym = create_sum(symbols)
        if new_sym not in const_dict:
            const_dict[new_sym] = len(const_dict)
        s = create_symbol(f_name(const_dict[new_sym]), sym_type)
        return [s]
    return symbols

