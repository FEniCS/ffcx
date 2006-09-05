__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2004-10-05 -- 2005-09-20"
__copyright__ = "Copyright (C) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

# FFC compiler modules
from referencetensor import *
from geometrytensor import *
import signature

class Term:
    """A Term represents a term of a multi-linear form and is
    expressed as the product of a ReferenceTensor A0 and a
    GeometryTensor GK.

    Attributes:

        product - the Product generating the Term
        A0      - reference tensor
        GKs     - list of geometry tensors (factorized)
    """

    def __init__(self, product, A0, GKs):
        "Create Term."

        self.product = product
        self.A0 = A0
        self.GKs = GKs

        return

    def signature(self):
        "Return signature of Term."
        return signature.compute_hard_signature(self.product)
