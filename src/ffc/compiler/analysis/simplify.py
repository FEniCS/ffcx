__author__ = "Marie Rognes (meg@math.uio.no)"
__date__ = "2006-10-23 -- 2007-03-20"
__copyright__ = "Copyright (C) 2006"
__license__  = "GNU GPL Version 2"

# Modified by Anders Logg 2007

# Python modules
import sys

# FFC common modules
from ffc.common.debug import *
from ffc.common.exceptions import *

# FFC compiler.language modules
from ffc.compiler.language.index import *
from ffc.compiler.language.algebra import *
from ffc.compiler.language.tokens import *

def simplify(form):
    """ Simplification of a Form f with respect to transforms and
    derivatives.

    This function will simplify terms reading:
    (dx_j/dX_i)(dX_l/dx_j) | (d/DX_l) => (d/dX_i)"""

    if not isinstance(form, Form):
        raise FormError, "I only know how to simplify a Form!"

    debug("Simplifying form...")

    # We aim to simplify each monomial on its own:
    for monomial in form.monomials:

        # Then we run thorough the basis functions in this monomial and
        # check for simplifications:
        for basis in monomial.basisfunctions:
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
                for transform in monomial.transforms:
                    if transform.type == Transform.JINV:
                        if not cmp(transform.index0, theindex):
                            first = transform
                            break
                for transform in monomial.transforms:
                    if transform.type == Transform.J:
                        if not cmp(transform.index1, first.index1):
                            second = transform
                            success = 1
                            break

                if success == 1:
                    # Now, we should first remove the transforms from
                    # the transform list. Second: replace the old
                    # derivative with the new one.
                    basis.derivatives[i] = Derivative(derivative.element, second.index0)
                    monomial.transforms.remove(first) 
                    monomial.transforms.remove(second)

    debug("done")
