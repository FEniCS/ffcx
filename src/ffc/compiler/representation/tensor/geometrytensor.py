__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2004-11-03 -- 2007-02-27"
__copyright__ = "Copyright (C) 2004-2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# FFC common modules
from ffc.common.debug import *

# FFC language modules
from ffc.compiler.language.index import *
from ffc.compiler.language.reassignment import *

# FFC tensor representation modules
from multiindex import *

class GeometryTensor:
    """This class represents the geometry tensor for a monomial term
    of a multilinear form"""

    def __init__(self, monomial):
        "Create geometry tensor for given monomial"

        debug("Computing geometry tensor...")

        # Save constants, coefficients and transforms
        self.constants = monomial.constants
        self.coefficients = monomial.coefficients
        self.transforms = monomial.transforms

        # Create secondary and auxiliary multi indices
        self.a = self.__create_multi_index(monomial, Index.SECONDARY)
        self.b = self.__create_multi_index(monomial, Index.AUXILIARY_G)
        debug("Secondary multi index: " + str(self.a), 1)
        debug("Auxiliary multi index: " + str(self.b), 1)

        debug("done")

    def __create_multi_index(self, monomial, index_type):
        "Find dimensions and create multi index"
        
        # Compute rank
        rank = max([max_index(c, index_type) for c in monomial.coefficients] + \
                   [max_index(t, index_type) for t in monomial.transforms] + [-1]) + 1
        
        # Compute all dimensions
        dims = [self.__find_dim(monomial, i, index_type) for i in range(rank)]

        # Create multi index from dims
        return MultiIndex(dims)

    def __find_dim(self, monomial, i, index_type):
        "Find dimension of given Index."

        # Create index to search for
        index = Index(i)
        index.type = index_type
        
        # Check coefficients
        for c in monomial.coefficients:
            if c.index == index:
                return c.e1.space_dimension()
            
        # Check transforms
        for t in monomial.transforms:
            if t.index0 == index or t.index1 == index:
                return t.element.cell_dimension()
            
        # Didn't find dimension
        raise RuntimeError, "Unable to find dimension for index " + str(index)
