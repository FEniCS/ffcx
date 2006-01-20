__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2004-11-06 -- 2005-11-08"
__copyright__ = "Copyright (c) 2004, 2005 Anders Logg"
__license__  = "GNU GPL Version 2"

# Python modules
from sets import Set # (Could use built-in set with Python 2.4)

# FFC common modules
from ffc.common.debug import *
from ffc.common.constants import *
from ffc.common.exceptions import *

# FFC compiler modules
from term import *
from reorder import *
from declaration import *

class ElementTensor:
    """An ElementTensor represents the element tensor of a
    multi-linear form and consist of a list of Terms, each containing
    a pair of a ReferenceTensor and a GeometryTensor.

    Attributes:

        terms - a list of Terms (products A0 * GK)
        a0    - a list of precomputed reference tensor declarations
        gK    - a list of precomputed geometry tensor declarations
        aK    - a list of precomputed element tensor declarations
    """

    def __init__(self, sum, type, format, cK_used):
        "Create ElementTensor."

        # Check that all Products have integrals
        self.__check_integrals(sum)

        # Reorder indices and compute factorization
        factorization = reorder_indices(sum)

        # Compute terms
        self.terms = [None for i in range(len(sum.products))]
        for i in range(len(sum.products)):
            debug("Compiling term %d" % i, 1)
            p = sum.products[i]
            if p.integral.type == type:
                # Compute geometry tensor
                GK = GeometryTensor(p)
                # Check if reference tensor should be computed
                if factorization[i] == None:
                    # Compute reference tensor and add term
                    A0 = ReferenceTensor(p)
                    self.terms[i] = Term(p, A0, [GK])
                else:
                    # Add geometry tensor to previous term
                    self.terms[factorization[i]].GKs += [GK]
        debug("All terms compiled", 1)

        # Remove terms not computed (factorized)
        [self.terms.remove(None) for i in range(len(self.terms)) if None in self.terms]
    
        # Compute reference tensor declarations (just need to pick the values)
        self.a0 = self.__compute_reference_tensor(format)

        # Compute element tensor declarations
        self.aK = self.__compute_element_tensor(format)

        # Compute geometry tensor declarations
        gK_used = self.__check_used(format)
        self.gK = self.__compute_geometry_tensor(format, gK_used, cK_used)

        return

    def __compute_reference_tensor(self, format):
        "Precompute reference tensor according to given format."
        debug("Generating code for reference tensor", 1)
        if not self.terms or format.format["reference tensor"](0, 0, []) == None: return []
        declarations = []
        for j in range(len(self.terms)):
            term = self.terms[j]
            iindices = term.A0.i.indices or [[]]
            aindices = term.A0.a.indices or [[]]
            for i in iindices:
                for a in aindices:
                    name = format.format["reference tensor"](j, i, a)
                    value = format.format["floating point"](term.A0(i, a))
                    declarations += [Declaration(name, value)]
        return declarations

    def __compute_geometry_tensor(self, format, gK_used, cK_used):
        "Precompute geometry tensor according to given format."
        debug("Generating code for geometry tensor", 1)
        if not self.terms or format.format["geometry tensor"](0, []) == None: return []
        declarations = []
        for j in range(len(self.terms)):
            # Should be the same, so pick first
            aindices = self.terms[j].GKs[0].a.indices
            if not aindices:
                aindices = [[]]
            for a in aindices:
                # Sum factorized values
                name = format.format["geometry tensor"](j, a)
                used = name in gK_used
                value = format.format["sum"]([GK(a, format, cK_used, used) for GK in self.terms[j].GKs])
                declaration = Declaration(name, value)
                declaration.used = used
                # Add declaration
                declarations += [declaration]

        return declarations

    def __compute_element_tensor(self, format):
        """Precompute element tensor, including optimizations. This is
        where any FErari optimization should be done."""
        debug("Generating code for element tensor", 1)
        if not self.terms or format.format["element tensor"]((0,), 0) == None: return []
        declarations = []
        iindices = self.terms[0].A0.i.indices or [[]] # All primary ranks are equal
        k = 0 # Update counter for each entry of A0, which is needed for some formats
        num_dropped = 0
        for i in iindices:
            debug("i = " + str(i), 2)
            value = None
            for j in range(len(self.terms)):
                debug("  j = " + str(j), 2)
                A0 = self.terms[j].A0
                if A0.a.indices: aindices = A0.a.indices
                else: aindices = [[]]
                for a in aindices:
                    debug("    a = " + str(a), 2)
                    name = format.format["element tensor"](i, k)
                    a0 = A0(i, a)
                    gk = format.format["geometry tensor"](j, a)
                    debug("      a0 = " + str(a0), 2)
                    debug("      gk = " + str(gk), 2)
                    if abs(a0) > FFC_EPSILON:
                        if value and a0 < 0.0:
                            value = format.format["subtract"]([value, format.format["multiplication"]([format.format["floating point"](-a0), gk])])
                        elif value:
                            value = format.format["sum"]([value, format.format["multiplication"]([format.format["floating point"](a0), gk])])
                        else:
                            value = format.format["multiplication"]([format.format["floating point"](a0), gk])
                    else:
                        num_dropped += 1
            value = value or format.format["floating point"](0.0)
            declarations += [Declaration(name, value)]
            k += 1
        debug("Number of zeros dropped from reference tensor: " + str(num_dropped), 1)
        return declarations

    def __check_used(self, format):
        """Check which declarations of gK are actually used, i.e,
        which entries of the geometry tensor that get multiplied with
        nonzero entries of the reference tensor."""
        gK_used = Set()
        if not self.terms or format.format["geometry tensor"](0, []) == None: return []
        iindices = self.terms[0].A0.i.indices or [[]] # All primary ranks are equal
        for i in iindices:
            for j in range(len(self.terms)):
                A0 = self.terms[j].A0
                if A0.a.indices: aindices = A0.a.indices
                else: aindices = [[]]
                for a in aindices:
                    a0 = A0(i, a)
                    gk = format.format["geometry tensor"](j, a)
                    if abs(a0) > FFC_EPSILON:
                        gK_used.add(gk)
        return gK_used

    def __check_integrals(self, sum):
        "Check that all terms have integrals."
        for p in sum.products:
            if not p.integral:
                raise FormError, (p, "Missing integral in term.")
