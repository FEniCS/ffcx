__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2004-11-03 -- 2007-03-05"
__copyright__ = "Copyright (C) 2004-2007 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Garth N. Wells 2006

# Python modules
import time
import numpy

# FFC language modules
#from ffc.compiler.language.index import *
#from ffc.compiler.language.reassignment import *

# FFC tensor representation modules
from uflmonomialintegration import integrate
#from multiindex import *

class ReferenceTensor:
    "This class represents the reference tensor for a monomial."

    def __init__(self, monomial, domain_type, facet0, facet1):
        "Create reference tensor for given monomial."

        # Compute reference tensor
        self.A0 = integrate(monomial, domain_type, facet0, facet1)

#        # Create primary, secondary and auxiliary multi indices
#        self.i = self.__create_multi_index(monomial, Index.PRIMARY)
#        self.a = self.__create_multi_index(monomial, Index.SECONDARY)
#        self.b = self.__create_multi_index(monomial, Index.AUXILIARY_0)
#        debug("Primary multi index: " + str(self.i), 1)
#        debug("Secondary multi index: " + str(self.a), 1)
#        debug("Auxiliary multi index: " + str(self.b), 1)

    def __create_multi_index(self, monomial, index_type):
        "Find dimensions and create multi index"
        
        # Compute rank
        rank = max([max_index(v, index_type) for v in monomial.basisfunctions] + [-1]) + 1

        # Compute all dimensions
        dims = [self.__find_dim(monomial, i, index_type) for i in range(rank)]

        # Create multi index from dims
        return MultiIndex(dims)

    def __find_dim(self, monomial, i, index_type):
        "Find dimension of given index"

        # Create index to search for
        index = Index(i)
        index.type = index_type

        # Search all basis functions
        for v in monomial.basisfunctions:
            
            # Check basis function index
            if v.index == index:
                return range(len(v.index.range))

            # Check component indices
            for j in range(len(v.component)):
                if v.component[j] == index:
                    return range(len(v.component[j].range))

            # Check derivatives
            for d in v.derivatives:
                if d.index == index:
                    return range(len(d.index.range))
                
        # Didn't find dimension
        raise RuntimeError, "Unable to find dimension for index " + str(index)
