__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2004-11-03 -- 2009-03-08"
__copyright__ = "Copyright (C) 2004-2009 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Marie E. Rognes (meg@math.uio.no) 2007

# FFC common modules
from ffc.common.log import debug

# FFC language modules
#from ffc.compiler.language.index import *
#from ffc.compiler.language.reassignment import *

# FFC tensor representation modules
from monomialextraction import MonomialException
from monomialtransformation import MonomialIndex
from multiindex import MultiIndex

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
        self.a = _create_multi_index(monomial, MonomialIndex.SECONDARY)
        self.b = _create_multi_index(monomial, MonomialIndex.AUXILIARY)

        debug("Secondary multi index: " + str(self.a))
        debug("Auxiliary multi index: " + str(self.b))

def _create_multi_index(monomial, index_type):
    "Find dimensions and create multi index for given index type."

    # Get sorted unique monomial indices
    indices = []
    for index in monomial.indices():
        if index.index_type == index_type and not index in indices:
            indices.append(index)
    indices.sort()

    # Check that we got all indices correctly
    for (i, index) in enumerate(indices):
        if not i == index.index_id:
            raise MonomialException, "Unable to extract all indices."

    # Get dimensions
    dims = [range(len(index.index_range)) for index in indices]

    return MultiIndex(dims)
