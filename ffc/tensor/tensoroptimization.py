__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2010-02-08"
__copyright__ = "Copyright (C) 2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# FFC modules
from ffc.log import warning

def optimize_integral_ir(ir):
    """
    Compute optimized intermediate representation of integral.

    Note that this function modifies the given intermediate
    representation directly, rather than working on a copy.
    """
    # Try importing FErari
    try:
        from ferari import binary
    except:
        warning("Unable to find FErari, skipping tensor optimizations")
        return ir

    # Iterate over tensors
    for (i, (A0, GK, optimized_contraction)) in enumerate(ir["AK"]):

        # Compute optimized tensor contraction
        if ir["rank"] == 2:
            optimized_contraction = binary.optimize(A0.A0)
        elif ir["rank"] == 1:
            optimized_contraction = binary.optimize_action(A0.A0)
        else:
            warning("Tensor optimization only available for rank 1 and 2 tensors, skipping optimizations")
            return ir

        # Store result for later use in code generation
        ir["AK"][i] = (A0, GK, optimized_contraction)

    return ir
