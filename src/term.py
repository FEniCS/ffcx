__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2004-10-13"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version "2

class Term:

    """A Term holds auxiliary data for a Product, including
    the rank and dimension of the Product.

    A Term holds the following data:

        r0   - primary tensor rank
        r1   - secondary tensor rank
        dims - a list of dimensions of length r0 + r1"""

    def __init__(self):
        self.r0 = 0
        self.r1 = 0
        dims = []
        return
