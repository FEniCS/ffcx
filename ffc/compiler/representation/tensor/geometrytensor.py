__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2004-11-03 -- 2009-03-08"
__copyright__ = "Copyright (C) 2004-2009 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Marie E. Rognes (meg@math.uio.no) 2007

# FFC common modules
from ffc.common.log import debug

# FFC tensor representation modules
from monomialtransformation import MonomialIndex
from multiindex import create_multi_index

class GeometryTensor:
    """This class represents the geometry tensor for a monomial term
    of a multilinear form."""

    def __init__(self, monomial):
        "Create geometry tensor for given monomial."

        # Save monomial data
        self.determinant = monomial.determinant
        self.coefficients = monomial.coefficients
        self.transforms = monomial.transforms

        # Create secondary and auxiliary multi indices
        self.secondary_multi_index = create_multi_index(monomial, MonomialIndex.SECONDARY)
        self.external_multi_index  = create_multi_index(monomial, MonomialIndex.EXTERNAL)

        debug("Secondary multi index: " + str(self.secondary_multi_index))
        debug("External multi index:  " + str(self.external_multi_index))
