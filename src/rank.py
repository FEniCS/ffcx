__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2004-10-13"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

# FFC modules
from index import *
from algebra import *
from reassign import *

class Rank:

    """A Rank holds auxiliary data for a Product, including
    the rank and dimension of the Product.

    A Rank holds the following data:

        r0    - primary tensor rank
        r1    - secondary tensor rank
        r2    - auxiliary tensor rank
        dims0 - list of dimensions for primary tensor indices
        dims1 - list of dimensions for secondary tensor indices
        dims2 - list of dimensions for auxiliary tensor indices"""
    
    def __init__(self, product):

        # Check that we get a Product
        if isinstance(product, Rank):
            self.r0 = product.r0
            self.r1 = product.r1
            self.r2 = product.r2
            self.dims0 = [] + product.dims0
            self.dims1 = [] + product.dims1
            self.dims2 = [] + product.dims2
            self.indices0 = [] + product.indices0
            self.indices1 = [] + product.indices1
            self.indices2 = [] + product.indices2
        elif isinstance(product, Product):

            # Compute ranks (works even if there is no index)
            self.r0 = max_product(product, "primary") + 1
            self.r1 = max_product(product, "secondary") + 1
            self.r2 = max_product(product, "auxiliary") + 1

            # Compute dimensions
            self.dims0 = dims_product(product, self.r0, "primary")
            self.dims1 = dims_product(product, self.r1, "secondary")
            self.dims2 = dims_product(product, self.r2, "auxiliary")

            # Compute index combinations
            self.indices0 = build_indices(self.dims0)
            self.indices1 = build_indices(self.dims1)
            self.indices2 = build_indices(self.dims2)

        else:
            raise RuntimeError, "Rank can only be created for a Product."

        print "Ranks: " + str(self.r0) + ", " + str(self.r1) + ", " + str(self.r2) + \
              " dimensions: " + str(self.dims0) + ", " + str(self.dims1) + ", " + str(self.dims2)
