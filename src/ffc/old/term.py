__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2004-10-05 -- 2007-01-23"
__copyright__ = "Copyright (C) 2004-2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# FFC formlang modules
from ffc.formlang.signature import *

# FFC compiler modules
from referencetensor import *
from geometrytensor import *

class Term:
    """A Term represents a term of a multi-linear form and is
    expressed as the product of a ReferenceTensor A0 and a
    GeometryTensor G.

    Attributes:

        monomial - the Monomial generating the Term
        A0       - reference tensor
        G        - list of geometry tensors (factorized)
    """

    def __init__(self, product, A0, G):
        "Create Term."

        self.product = product
        self.A0 = A0
        self.G  = G

        return

    def signature(self):
        "Return signature of Term."
        return compute_hard_signature(self.product)
