__author__ = "Kristian B. Oelgaard (k.b.oelgaard@gmail.com)"
__date__ = "2010-02-08"
__copyright__ = "Copyright (C) 2010 Kristian B. Oelgaard"
__license__  = "GNU GPL version 3 or any later version"

# FFC modules
from ffc.log import info

def optimize_integral_ir(ir):
    "Compute optimized intermediate representation of integral."

    info("Optimisation of quadrature code takes place at the code generation stage, skipping optimization")

    return ir
