__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2004-11-03"
__copyright__ = "Copyright (C) 2004-2009 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Marie E. Rognes (meg@math.uio.no) 2007
# Modified by Kristian B. Oelgaard, 2009
# Last changed: 2009-12-21

# FFC modules.
from ffc.log import debug

# FFC tensor representation modules.
from monomialtransformation import MonomialIndex
from multiindex import create_multiindex

class GeometryTensor:
    """
    This class represents the geometry tensor for a monomial term of a
    multilinear form.
    """

    def __init__(self, monomial):
        "Create geometry tensor for given monomial."

        # Save monomial data
        self.determinant = monomial.determinant
        self.coefficients = monomial.coefficients
        self.transforms = monomial.transforms

        # Extract indices
        secondary_indices = monomial.extract_unique_indices(MonomialIndex.SECONDARY)
        external_indices  = monomial.extract_unique_indices(MonomialIndex.EXTERNAL)

        # Create multiindices
        self.secondary_multi_index = create_multiindex(secondary_indices)
        self.external_multi_index  = create_multiindex(external_indices)

        debug("Secondary multi index: " + str(self.secondary_multi_index))
        debug("External multi index:  " + str(self.external_multi_index))
