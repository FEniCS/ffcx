__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2004-10-05 -- 2005-09-07"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

# FFC modules
from referencetensor import *
from geometrytensor import *

class Term:

    """A Term represents a term of a multi-linear form and is
    expressed as the product of a ReferenceTensor A0 and a
    GeometryTensor GK.

        A0  - reference tensor
        GKs - list of geometry tensors (factorized)"""

    def __init__(self, A0, GKs):
        "Create Term."

        self.A0 = A0
        self.GKs = GKs

        return
