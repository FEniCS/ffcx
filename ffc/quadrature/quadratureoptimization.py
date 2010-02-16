__author__ = "Kristian B. Oelgaard (k.b.oelgaard@gmail.com)"
__date__ = "2010-02-08"
__copyright__ = "Copyright (C) 2010 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# FFC modules
from ffc.log import info
from ffc.quadrature.symbolics import optimise_code

def optimize_integral_ir(ir):
    "Compute optimized intermediate representation of integral."

    parameters = ir["optimise_parameters"]
    if parameters["simplify expressions"]:
        integrals  = ir["trans_integrals"]
        domain_type =ir["domain_type"]
        num_facets = ir["num_facets"]
        geo_consts = ir["geo_consts"]
        if domain_type == "cell":
            info("Optimising expressions for cell integral")
            _optimise_integral(integrals, geo_consts)
        elif domain_type == "exterior_facet":
            for i in range(num_facets):
                info("Optimising expressions for facet integral %d" % i)
                _optimise_integral(integrals[i], geo_consts)
        elif domain_type == "interior_facet":
            for i in range(num_facets):
                for j in range(num_facets):
                    info("Optimising expressions for facet integral (%d, %d)" % (i, j))
                    _optimise_integral(integrals[i][j], geo_consts)
        else:
            error("Unhandled domain type: " + str(domain_type))

    return ir

def _optimise_integral(integral, geo_consts):
    for points, terms, functions, ip_consts in integral:
        new_terms = {}
        for key, data in terms.iteritems():
            val, ops, t_set, u_weights, u_psi_tables, u_nzcs = data
            value = optimise_code(val, ip_consts, geo_consts, t_set)
            # Check if value is zero
            if value.val:
                data[0] = value
                data[1] = value.ops()
                data[2] = t_set
                new_terms[key] = data
        terms = new_terms

