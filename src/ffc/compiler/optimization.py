__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2006-03-22"
__copyright__ = "Copyright (C) 2006 Anders Logg"
__license__  = "GNU GPL Version 2"

# FFC common modules
from ffc.common.debug import *

def optimize(A0):
    """Generate optimized abstract code for the tensor contraction from
    a given reference tensor"""

    debug("Computing optimization, this may take some time...")

    # Check if FErari is available
    try:
        from ferari import binary
    except:
        raise RuntimeError, "Cannot find FErari on your system, unable to optimize."

    print ""
    print "--- Reference tensor ---"
    print A0
    print ""

    # Compute optimized code
    code = binary.optimize(A0)

    print "--- Optimized code ---"
    print code
    print ""
