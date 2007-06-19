"Element tensor class for quadrature representation"

__author__ = "Kristian B. Oelgaard (k.b.oelgaard@tudelft.nl)"
__date__ = "2007-03-16 -- 2007-06-19"
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

        monomial               - monomial from which the element was constructed
        Derivatives            - table containing derivatives (to construct Jacobian)
        Psis                   - table containing the values of the basisfunctions
        quadrature             - Quadrature object, (points, weights)
        i,a,b0,bg              - indices
        constants
        coefficients
        coefficient_offsets    - offsets for restricted elements
        transforms
        determinant
    """

    def __init__(self, monomial, facet0, facet1):
        "Create element tensor for given monomial"

        # Save monomial
        self.monomial = monomial

        # Tabulate element tensor
        self.Psis, self.quadrature = tabulate(monomial, facet0, facet1)

        # Create primary, secondary and auxiliary multi indices
        self.i = self.__create_multi_index(monomial, Index.PRIMARY)
        self.a = self.__create_multi_index(monomial, Index.SECONDARY)
        self.b0 = self.__create_multi_index(monomial, Index.AUXILIARY_0)
        self.bg = self.__create_multi_index(monomial, Index.AUXILIARY_G)

        # Save constants, coefficients and transforms (used to generate factor)
        self.constants = monomial.constants
        self.coefficients = monomial.coefficients

        # Compute offsets (for restricted elements)
        self.coefficient_offsets = self.__compute_coefficient_offsets(monomial)
        self.transforms = monomial.transforms
        self.determinant = monomial.determinant

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

        # Initialise dictionary
        offsets = {}

        # Compute offset dependent on restriction and range of coefficient
        for c in monomial.coefficients:
            for v in monomial.basisfunctions:
                if c.index == v.index:
                    offsets[c] = choose_map[v.restriction]*len(c.index.range)
                    break
        return offsets


