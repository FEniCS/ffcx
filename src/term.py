__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2004-10-05"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

# FFC modules
from referencetensor import *
from geometrytensor import *

class Term:

    """A Term represents a term of a multi-linear form and is
    expressed as the product of a ReferenceTensor A0 and a
    GeometryTensor GK."""

    def __init__(self, product):
        "Create Term."
        # Create ReferenceTensor and GeometryTensor
        self.A0 = ReferenceTensor(product)
        self.GK = GeometryTensor(product)
        # Check that the ranks match
        if not self.A0.a.rank == self.GK.a.rank:
            raise RuntimeError, "Secondary ranks don't match."            
        return
