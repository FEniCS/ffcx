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
                return t.element.shapedim()
            
        # Didn't find dimension
        raise RuntimeError, "Unable to find dimension for index " + str(index)

# from util import *

# # FFC compiler.language modules
# from algebra import *
# from reassign import *
# from multiindex import *

# class GeometryTensor:

#     """A GeometryTensor represents the geometry tensor of a
#     multi-linear form computed on the current cell.

#         a            - secondary multiindex
#         b            - auxiliary multiindex
#         rank         - rank of the tensor
#         constants    - a list of Constants
#         coefficients - a list of Coefficients
#         transforms   - a list of Transforms"""

#     def __init__(self, monomial):
#         "Create GeometryTensor."

#         # Check that we get a Monomial
#         if not isinstance(monomial, Monomial):
#             raise RuntimeError, "GeometryTensor must be created from Monomial."

#         # Get data from Monomial
#         self.constants = listcopy(monomial.constants)
#         self.coefficients = listcopy(monomial.coefficients)
#         self.transforms = listcopy(monomial.transforms)

#         # Get rank
#         self.rank = self.a.rank

#         debug("Created geometry tensor: a%s, b%s" % \
#               (str(self.a.dims), str(self.b.dims)), 1)

#         return


#     def __call__(self, a, format, cK_used, used):
#         "Return given element of geometry tensor."
#         # Compute product of factors outside sum
#         factors = []
#         for c in self.constants:
#             if c.inverted:
#                 factors += ["(1.0/" + format.format["constant"](c.number.index) + ")"]
#             else:
#                 factors += [format.format["constant"](c.number.index)]
#         for c in self.coefficients:
#             if not c.index.type == Index.AUXILIARY_G:
#                 coefficient = format.format["coefficient"](c.n1.index, c.index([], a, [], []))
#                 factors += [coefficient]
#                 # Only add coefficients appearing in an entry of G that is used
#                 if used: cK_used.add(coefficient)
#         for t in self.transforms:
#             if not (t.index0.type == Index.AUXILIARY_G or  t.index1.type == Index.AUXILIARY_G):
#                 factors += [format.format["transform"](t.index0([], a, [], []), \
#                                                        t.index1([], a, [], []), \
#                                                        t.restriction),]
#         monomial = format.format["multiplication"](factors)
#         if monomial: f0 = [monomial]
#         else: f0 = []
#         # Compute sum of monomials inside sum
#         terms = []
#         for b in self.b.indices:
#             factors = []
#             for c in self.coefficients:
#                 if c.index.type == Index.AUXILIARY_G:
#                     coefficient = format.format["coefficient"](c.n1.index, c.index([], a, [], b))
#                     factors += [coefficient]
#                     # Only add coefficients appearing in an entry of G that is used
#                     if used: cK_used.add(coefficient)
#             for t in self.transforms:
#                 if t.index0.type == Index.AUXILIARY_G or t.index1.type == Index.AUXILIARY_G:
#                     factors += [format.format["transform"](t.index0([], a, [], b), \
#                                                            t.index1([], a, [], b), \
#                                                            t.restriction)]
#             terms += [format.format["multiplication"](factors)]
#         sum = format.format["sum"](terms)
#         if sum: sum = format.format["grouping"](sum)
#         if sum: f1 = [sum]
#         else: f1 = []
#         # Compute product of all factors
#         return format.format["multiplication"]([f for f in [format.format["determinant"]] + f0 + f1])

#     def __str__(self):
#         "Print nicely formatted representation of GeometryTensor."
#         c = "".join([c.__str__() for c in self.coefficients])
#         t = "".join([t.__str__() for t in self.transforms])
#         return c + t
