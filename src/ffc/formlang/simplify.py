__author__ = "Marie Rognes (meg@math.uio.no)"
__date__ = "2006-10-23 -- 2007-01-23"
__copyright__ = "Copyright (C) 2006"
__license__  = "GNU GPL Version 2"

# Modified by Anders Logg 2007

# Python modules
import sys

# FFC common modules
sys.path.append("../")
from ffc.common.debug import *
from ffc.common.exceptions import *

# FFC formlang modules
from index import *
from algebra import *
from tokens import *

def simplify(f):
    """ Simplification of a Sum f with respect to transforms and
    derivatives.

    This function will simplify terms reading:
    (dx_j/dX_i)(dX_l/dx_j) | (d/DX_l) => (d/dX_i)"""

    if not isinstance(f, Sum):
        raise FormError, "I only know how to simplify a Sum!"

    modified = Sum(f)

    # We aim to simplify each product on its own:
    for product in modified.products:

        # Then we run thorough the basis functions in this product and
        # check for simplifications:
        for basis in product.basisfunctions:
            success = 0
            first = None
            second = None

            # The derivatives of this basis function is a good
            # starting point so we run through these. We use index
            # notation since we may have to replace some of them.
            for i in range(len(basis.derivatives)):
                derivative = basis.derivatives[i]
                success = 0
                theindex = derivative.index
                # Now, lets run through the transforms and see whether
                # there are two matching:
                for transform in product.transforms:
                    # For simplicity assume that transform.power == 1 | -1
                    if transform.power == 1:
                        if not cmp(transform.index0, theindex):
                            first = transform
                            break
                for transform in product.transforms:
                    if transform.power == -1:
                        if not cmp(transform.index1, first.index1):
                            second = transform
                            success = 1
                            break

                if success == 1:
                    # Now, we should first remove the transforms from
                    # the transform list. Second: replace the old
                    # derivative with the new one.
                    basis.derivatives[i] = Derivative(derivative.element, second.index0)
                    product.transforms.remove(first) 
                    product.transforms.remove(second)

    return modified
