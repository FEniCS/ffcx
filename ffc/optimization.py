"""
Compiler stage 5: optimization
------------------------------

This module implements the optimization of an intermediate code
representation.
"""

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2009-12-22"
__copyright__ = "Copyright (C) 2009 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Last changed: 2010-02-08

# FFC modules
from ffc.log import info, begin, end

# FFC specialized code generation modules
from ffc import quadrature
from ffc import tensor

def optimize_ir(ir, parameters):
    "Optimize intermediate form representation."

    begin("Compiler stage 3: Optimizing intermediate representation")

    # Check if optimization is requested
    if not parameters["optimize"]:
        info("Skipping optimizations, add -O to optimize")
        return ir

    # Extract representations
    ir_elements, ir_dofmaps, ir_integrals, ir_forms = ir

    # Iterate over integrals
    oir_integrals = [_optimize_integral_ir(ir) for ir in ir_integrals]

    return ir_elements, ir_dofmaps, oir_integrals, ir_forms

def _optimize_integral_ir(ir):
    "Compute optimized intermediate represention of integral."

    # Select representation
    if ir["representation"] == "quadrature":
        r = quadrature
    elif ir["representation"] == "tensor":
        r = tensor
    else:
        error("Unknown representation: %s" % ir["representation"])

    # Optimize representation
    oir = r.optimize_integral_ir(ir)

    return ir
