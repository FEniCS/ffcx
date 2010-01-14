__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2004-11-03"
__copyright__ = "Copyright (C) 2004-2009 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Garth N. Wells 2006
# Modified by Kristian B. Oelgaard, 2009.
# Last changed: 2010-01-14

# FFC modules.
from ffc.log import debug

# FFC tensor representation modules.
from monomialintegration import integrate
from monomialtransformation import MonomialIndex
from multiindex import create_multiindex

class ReferenceTensor:
    """
    This class represents the reference tensor for a monomial term of
    a multilinear form.
    """

    def __init__(self, monomial, domain_type, facet0, facet1, quadrature_order):
        "Create reference tensor for given monomial."

        # Compute reference tensor
        self.A0 = integrate(monomial, domain_type, facet0, facet1, quadrature_order)

        # Extract indices
        primary_indices   = monomial.extract_unique_indices(MonomialIndex.PRIMARY)
        secondary_indices = monomial.extract_unique_indices(MonomialIndex.SECONDARY)
        internal_indices  = monomial.extract_unique_indices(MonomialIndex.INTERNAL)

        # Create multiindices
        self.primary_multi_index   = create_multiindex(primary_indices)
        self.secondary_multi_index = create_multiindex(secondary_indices)
        self.internal_multi_index  = create_multiindex(internal_indices)

        # Store monomial
        self.monomial = monomial

        debug("Primary multi index:   " + str(self.primary_multi_index))
        debug("Secondary multi index: " + str(self.secondary_multi_index))
        debug("Internal multi index:  " + str(self.internal_multi_index))
