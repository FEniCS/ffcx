__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2004-10-13"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

from algebra import *
from reassign import *

class Rank:

    """A Rank holds auxiliary data for a Product, including
    the rank and dimension of the Product.

    A Rank holds the following data:

        r0   - primary tensor rank
        r1   - secondary tensor rank
        dims - a list of dimensions of length r0 + r1"""

    def __init__(self, product):

        # Check that we get a Product
        if isinstance(product, Rank):
            self.r0 = product.r0
            self.r1 = product.r1
            self.dims = [] + product.dims
        elif isinstance(product, Product):
            # Note that this works even if there is no index.
            # In that case, we get -1 + 1 = 0, which is correct.
            self.r0 = max_product(product, "primary") + 1
            self.r1 = max_product(product, "secondary") + 1
            self.dims = dims_product(product, self.r0, self.r1)
        else:
            raise RuntimeError, "Rank can only be created for a Product."

        print "Primary rank: " + str(self.r0) + \
              ", secondary rank: " + str(self.r1) + \
              ", dimensions: " + str(self.dims)
