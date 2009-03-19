__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2004-11-03 -- 2009-03-09"
__copyright__ = "Copyright (C) 2004-2009 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Garth N. Wells 2006

# FFC common modules
from ffc.common.log import debug

# FFC tensor representation modules
from uflmonomialintegration import integrate
from monomialtransformation import MonomialIndex
from multiindex import create_multi_index

class ReferenceTensor:
    "This class represents the reference tensor for a monomial."

    def __init__(self, monomial, domain_type, facet0, facet1):
        "Create reference tensor for given monomial."

        # Compute reference tensor
        self.A0 = integrate(monomial, domain_type, facet0, facet1)

        # FIXME: Handle auxiliary indices for A0 and GK

        # Create primary, secondary and auxiliary multi indices
        self.i = create_multi_index(monomial, MonomialIndex.PRIMARY)
        self.a = create_multi_index(monomial, MonomialIndex.SECONDARY)
        self.b = create_multi_index(monomial, MonomialIndex.AUXILIARY)

        debug("Primary multi index: " + str(self.i))
        debug("Secondary multi index: " + str(self.a))
        debug("Auxiliary multi index: " + str(self.b))

        import sys
        sys.exit(1)
