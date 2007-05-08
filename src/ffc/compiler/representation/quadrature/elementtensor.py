__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2007-03-16 -- 2007-03-23"
__copyright__ = "Copyright (C) 2007 Kristian B. Oelgaard"
__license__  = "GNU GPL Version 2"

# Python modules
import time
import numpy

# FFC language modules
from ffc.compiler.language.index import *
from ffc.compiler.language.reassignment import *

# FFC tensor representation modules
#from monomialintegration import *
from monomialtabulation import *

from multiindex import *

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

        # Tabulate element tensor
        self.map_derivatives, self.map_element, self.Psis, self.quadrature = \
        tabulate(monomial, facet0, facet1)

        # FIXME: needed?
        # Create primary, secondary and auxiliary multi indices
        self.i = self.__create_multi_index(monomial, Index.PRIMARY)
#        print "element tensor, self.i : ", self.i
        self.a = self.__create_multi_index(monomial, Index.SECONDARY)
#        print "element tensor, self.a : ", self.a
        self.b = self.__create_multi_index(monomial, Index.AUXILIARY_0)
#        print "element tensor, self.b : ", self.b

        # Save constants, coefficients and transforms (used to generate factor)
        self.constants = monomial.constants
#        print "geo: constants: ", monomial.constants

        self.coefficients = monomial.coefficients
#        print "geo: coeff: ", monomial.coefficients
        self.transforms = monomial.transforms
#        print "geo: transforms: ", monomial.transforms

        self.determinant = monomial.determinant
#        print "geo: determinant: ", monomial.determinant

#        debug("Primary multi index: " + str(self.i), 1)
#        debug("Secondary multi index: " + str(self.a), 1)
#        debug("Auxiliary multi index: " + str(self.b), 1)

    # FIXME: needed?
    def __create_multi_index(self, monomial, index_type):
        "Find dimensions and create multi index"
        
        # Compute rank
        rank = max([max_index(v, index_type) for v in monomial.basisfunctions] + [-1]) + 1

        # Compute all dimensions
        dims = [self.__find_dim(monomial, i, index_type) for i in range(rank)]

        # Create multi index from dims
        return MultiIndex(dims)

    # FIXME: needed?
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
