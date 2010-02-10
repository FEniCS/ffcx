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

    warning("Calling FErari to perform optimizations, currently broken")

    # Iterate over tensors
    for (A0, GK) in ir["AK"]:

        # Compute optimized abstract code
        if ir["rank"] == 2:
            abstract_code = binary.optimize(A0.A0)
        elif ir["rank"] == 1:
            abstract_code = binary.optimize_action(A0.A0)
        else:
            warning("Tensor optimization only available for rank 1 and 2 tensors, skipping optimizations")
            return ir

        print abstract_code

    return ir
