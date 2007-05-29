__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2007-03-16 -- 2007-03-23"
__copyright__ = "Copyright (C) 2007 Kristian B. Oelgaard"
__license__  = "GNU GPL Version 2"

# FFC language modules
from ffc.compiler.language.index import *
from ffc.compiler.language.reassignment import *
from ffc.compiler.language.restriction import *

# FFC tensor representation modules
from ffc.compiler.representation.tensor.multiindex import *

from monomialtabulation import *

# Choose map from restriction, (from ufcformat.py - not sure about CONSTANT)
choose_map = {Restriction.PLUS: 0, Restriction.MINUS: 1, Restriction.CONSTANT: 0, None: 0}

class ElementTensor:
    """This class represents the element tensor for a monomial term
    of a multilinear form.

    Attributes:

        Derivatives            - table containing derivatives (to construct Jacobian)
        Psis                   - table containing the values of the basisfunctions
        quadrature             - Quadrature object, (points, weights)
    """

    def __init__(self, monomial, facet0, facet1):
        "Create element tensor for given monomial"

        # Save monomial
        self.monomial = monomial

        # Tabulate element tensor
        self.Psis, self.quadrature = tabulate(monomial, facet0, facet1)

#        print "element tensor, monomial : ", monomial
        # Create primary, secondary and auxiliary multi indices
        self.i = self.__create_multi_index(monomial, Index.PRIMARY)
#        print "element tensor, self.i : ", self.i
        self.a = self.__create_multi_index(monomial, Index.SECONDARY)
#        print "element tensor, self.a : ", self.a
        self.b0 = self.__create_multi_index(monomial, Index.AUXILIARY_0)
#        print "b0", self.b0
        self.bg = self.__create_multi_index(monomial, Index.AUXILIARY_G)
#        print "bg", self.bg

        # Save constants, coefficients and transforms (used to generate factor)
        self.constants = monomial.constants
#        print "geo: constants: ", monomial.constants

        self.coefficients = monomial.coefficients
#        print "geo: coeff: ", monomial.coefficients
        self.coefficient_offsets = self.__compute_coefficient_offsets(monomial)

        self.transforms = monomial.transforms
#        print "geo: transforms: ", monomial.transforms

        self.determinant = monomial.determinant
#        print "geo: determinant: ", monomial.determinant

        debug("Primary multi index: " + str(self.i), 1)
        debug("Secondary multi index: " + str(self.a), 1)
        debug("Auxiliary_0 multi index: " + str(self.b0), 1)
        debug("Auxiliary_G multi index: " + str(self.bg), 1)

    def __create_multi_index(self, monomial, index_type):
        "Find dimensions and create multi index"
        
        rank = max([max_index(v, index_type) for v in monomial.basisfunctions] + [-1]) + 1

        # Compute all dimensions (for reference tensor)
        dims = [self.__find_dim_ref(monomial, i, index_type) for i in range(rank)]

        # We're looking for Index.AUXILIARY_G
        if not dims:
            # Compute rank
            rank = max([max_index(c, index_type) for c in monomial.coefficients] + \
                       [max_index(t, index_type) for t in monomial.transforms] + [-1]) + 1

            # Compute all dimensions (for geometry tensor)
            dims = [self.__find_dim_geo(monomial, i, index_type) for i in range(rank)]

        # Create multi index from dims
        return MultiIndex(dims)

    def __find_dim_ref(self, monomial, i, index_type):
        "Find dimension of given index"

        # Create index to search for
        index = Index(i)
        index.type = index_type

        # Search all basis functions
        for v in monomial.basisfunctions:
            
            # Check basis function index
            if v.index == index:
#                return v.element.space_dimension()
                return range(len(v.index.range))

            # Check component indices
            for j in range(len(v.component)):
                if v.component[j] == index:
#                    return v.element.value_dimension(j)
                    return range(len(v.component[j].range))

            # Check derivatives
            for d in v.derivatives:
                if d.index == index:
#                    return d.element.cell_dimension()
                    return range(len(d.index.range))

                
        # Didn't find dimension
        raise RuntimeError, "Unable to find dimension for index " + str(index)

    def __find_dim_geo(self, monomial, i, index_type):
        "Find dimension of given index"

        # Create index to search for
        index = Index(i)
        index.type = index_type

        # Check coefficients
        for c in monomial.coefficients:
            if c.index == index:
                return range(len(c.index.range))
            
        # Check transforms
        for t in monomial.transforms:
            if t.index0 == index:
                return range(len(t.index0.range))
            elif t.index1 == index:
                return range(len(t.index1.range))

        # Didn't find dimension
        raise RuntimeError, "Unable to find dimension for index " + str(index)

    def __compute_coefficient_offsets(self, monomial):

        offsets = {}
#        print "mon: ", monomial
#        print "coeff: ", monomial.coefficients
#        for v in monomial.basisfunctions:
#            print "v: ", v
#            print "v.index: ", v.index
#            print "v.res: ", v.restriction

        for c in monomial.coefficients:
            for v in monomial.basisfunctions:
                if c.index == v.index:
                    offsets[c] = choose_map[v.restriction]*len(c.index.range)
#                    print "c: ", c
#                    print "c.index.range: ", c.index.range
#            print "c.n0: ", c.n0
#            print "c.n1: ", c.n1
#            print "c.e0: ", c.e0
#            print "c.e1: ", c.e1
#            print "c.P: ", c.P
#                    print "c.index: ", c.index
                    break
#            print "c.ops: ", c.ops
#        print offsets

#        raise RuntimeError, "restrict "
        return offsets


