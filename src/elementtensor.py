__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2004-11-06"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

EPSILON = 3e-16

# FFC modules
from term import *

class ElementTensor:

    """An ElementTensor represents the element tensor of a
    multi-linear form and consist of a list of Terms, each containing
    a pair of a ReferenceTensor and a GeometryTensor."""

    def __init__(self, sum, type, format):
        "Create ElementTensor."

        # Create list of Terms
        self.terms = [Term(p) for p in sum.products if p.integral.type == type]

        # Compute element tensor
        self.AK = self.__compute_element_tensor(format)
        
        return

    def __compute_element_tensor(self, format):
        """Precompute element tensor, including optimizations. This is
        where any FErari optimization should be done."""

        # Iterate over all primary indices and compute the values
        iindices = self.terms[0].A0.i.indices # All primary ranks are equal
        values = []
        for i in iindices:
            value = ""
            for j in range(len(self.terms)):
                A0 = self.terms[j].A0
                GK = self.terms[j].GK
                if A0.a.indices: aindices = A0.a.indices
                else: aindices = [[]]
                for a in aindices:
                    a0 = A0(i, a)
                    gk = GK.name(j, a, format)
                    if abs(a0) > EPSILON:
                        if value and a0 < 0.0:
                            value += " - %s*%s" % (str(-a0), gk)
                        elif value:
                            value += " + %s*%s" % (str(a0), gk)
                        else:
                            value += "%s*%s" % (str(a0), gk)
            values += [value]
            
        # Create dictionary with index and value pairs
        keys = [str(i) for i in iindices]
        return dict([(keys[j], values[j]) for j in range(len(keys))])

    def __call__(self, i):
        "Return given element of element tensor."
        return self.AK[str(i)]
