__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2004-11-06"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

EPSILON = 3e-16

# FFC modules
from term import *
from declaration import *

class ElementTensor:

    """An ElementTensor represents the element tensor of a
    multi-linear form and consist of a list of Terms, each containing
    a pair of a ReferenceTensor and a GeometryTensor.

    An ElementTensor holds the following data:

        terms  - a list of Terms (products A0 * GK)
        gK     - a list of precomputed geometry tensor declarations
        aK     - a list of precomputed element tensor declarations"""

    def __init__(self, sum, type, format):
        "Create ElementTensor."

        # Create list of Terms
        self.terms = [Term(p) for p in sum.products if p.integral.type == type]

        # Compute geometry tensor
        self.gK = self.__compute_geometry_tensor(format)

        # Compute element tensor
        self.aK = self.__compute_element_tensor(format)
        
        return

    def __compute_geometry_tensor(self, format):
        "Precompute geometry tensor according to given format."
        declarations = []
        for j in range(len(self.terms)):
            GK = self.terms[j].GK
            if GK.a.indices: aindices = GK.a.indices
            else: aindices = [[]]
            for a in aindices:
                name = format["geometry tensor"](j, a)
                value = GK(a, format)
                declarations += [Declaration(name, value)]
        return declarations

    def __compute_element_tensor(self, format):
        """Precompute element tensor, including optimizations. This is
        where any FErari optimization should be done."""
        aformat = format["element tensor"]
        iindices = self.terms[0].A0.i.indices # All primary ranks are equal
        declarations = []
        for i in iindices:
            value = ""
            for j in range(len(self.terms)):
                A0 = self.terms[j].A0
                GK = self.terms[j].GK
                if A0.a.indices: aindices = A0.a.indices
                else: aindices = [[]]
                for a in aindices:
                    name = format["element tensor"](i)
                    a0 = A0(i, a)
                    gk = format["geometry tensor"](j, a)
                    if abs(a0) > EPSILON:
                        if value and a0 < 0.0:
                            value += " - %s*%s" % (str(-a0), gk)
                        elif value:
                            value += " + %s*%s" % (str(a0), gk)
                        else:
                            value += "%s*%s" % (str(a0), gk)
            declarations += [Declaration(name, value)]
        return declarations    
