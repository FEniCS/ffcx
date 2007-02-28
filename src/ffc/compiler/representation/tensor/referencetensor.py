__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2004-11-03 -- 2007-02-27"
__copyright__ = "Copyright (C) 2004-2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# Modified by Garth N. Wells 2006

# Python modules
import time
import numpy

# FFC language modules
from ffc.compiler.language.index import *
from ffc.compiler.language.reassignment import *

# FFC tensor representation modules
from monomialintegration import *
from multiindex import *

class ReferenceTensor:
    """This class represents the reference tensor for a monomial term
    of a multilinear form"""

    def __init__(self, monomial, facet0, facet1):
        "Create reference tensor for given monomial"

        # Compute reference tensor
        self.A0 = integrate(monomial, facet0, facet1)

        # Create primary, secondary and auxiliary multi indices
        self.i = self.__create_multi_index(monomial, Index.PRIMARY)
        self.a = self.__create_multi_index(monomial, Index.SECONDARY)
        self.b = self.__create_multi_index(monomial, Index.AUXILIARY_0)
        debug("Primary multi index: " + str(self.i), 1)
        debug("Secondary multi index: " + str(self.a), 1)
        debug("Auxiliary multi index: " + str(self.b), 1)

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
                return v.element.space_dimension()

            # Check component indices
            for j in range(len(v.component)):
                if v.component[j] == index:
                    return v.element.value_dimension(j)

            # Check derivatives
            for d in v.derivatives:
                if d.index == index:
                    return d.element.cell_dimension()
                
        # Didn't find dimension
        raise RuntimeError, "Unable to find dimension for index " + str(index)

# # FFC common modules
# from exceptions import *
# from progress import *
# from debug import *
# from util import *

# # FFC compiler.language modules
# from index import *
# from algebra import *
# from reassign import *
# from multiindex import *



# class ReferenceTensor:

#     """A ReferenceTensor represents the reference tensor of a
#     multi-linear form computed on the reference cell.

#     A ReferenceTensor holds the following data:

#         i              - primary multiindex
#         a              - secondary multiindex
#         b              - auxiliary multiindex
#         rank           - rank of the tensor
#         A0             - the precomputed reference tensor
#         numeric        - a numeric constant (float)
#         basisfunctions - a list of BasisFunctions
#         integral       - an Integral
#         cputime        - time to compute the reference tensor"""

#     def __init__(self, monomial, facet0, facet1, alignment):
#         "Create ReferenceTensor."

#         # Check that we get a Monomial
#         if not isinstance(monomial, Monomial):
#             raise RuntimError, "ReferenceTensor must be created from Monomial."

#         # Check that the Monomial contains an Integral
#         if monomial.integral == None:
#             raise FormError, (monomial, "Missing integral in term.")

#         # Get data from Monomial
#         self.numeric = monomial.numeric
#         self.basisfunctions = listcopy(monomial.basisfunctions)
#         self.integral = monomial.integral



#         # Report time to compute the reference tensor
#         self.cputime = time.time() - t
#         debug("Reference tensor computed in %.3g seconds" % self.cputime)

#         # Get rank
#         self.rank = self.i.rank + self.a.rank

#         debug("Created reference tensor: i%s, a%s, b%s" % \
#               (str(self.i.dims), str(self.a.dims), str(self.b.dims)), 1)
        
#         return



#     def __find_dim(self, i, type):
#         "Find dimension of given Index."
#         index = Index(i)
#         index.type = type
#         for v in self.basisfunctions:
#             # Check basis Index
#             if v.index == index:
#                 return v.element.space_dimension()
#             # Check component Indices
#             for j in range(len(v.component)):
#                 if v.component[j] == index:
#                     return v.element.value_dimension(j)
#             # Check Derivatives
#             for d in v.derivatives:
#                 if d.index == index:
#                     return d.element.cell_dimension()
#         # Didn't find dimension
#         raise RuntimeError, "Unable to find dimension for Index " + str(index)

#     def  __call__(self, i, a = []):
#         "Return given element of reference tensor."
#         return self.A0[i + a]

#     def __str__(self):
#         "Pretty print"
#         v = "*".join([v.__str__() for v in self.basisfunctions])
#         return v + "*" + self.integral.__str__()
