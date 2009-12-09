__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2004-11-03"
__copyright__ = "Copyright (C) 2004-2009 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Garth N. Wells 2006
# Modified by Kristian B. Oelgaard, 2009.
# Last changed: 2009-12-09

# FFC modules.
from ffc.log import debug

# FFC tensor representation modules.
from monomialintegration import integrate
from monomialtransformation import MonomialIndex
from multiindex import create_multi_index

class ReferenceTensor:
    "This class represents the reference tensor for a monomial."

    def __init__(self, monomial, domain_type, facet0, facet1, quadrature_order):
        "Create reference tensor for given monomial."

        # Compute reference tensor
        self.A0 = integrate(monomial, domain_type, facet0, facet1, quadrature_order)

        # Create primary, secondary and auxiliary multi indices
        self.primary_multi_index   = create_multi_index(monomial, MonomialIndex.PRIMARY)
        self.secondary_multi_index = create_multi_index(monomial, MonomialIndex.SECONDARY)
        self.internal_multi_index  = create_multi_index(monomial, MonomialIndex.INTERNAL)

        debug("Primary multi index:   " + str(self.primary_multi_index))
        debug("Secondary multi index: " + str(self.secondary_multi_index))
        debug("Internal multi index:  " + str(self.internal_multi_index))
