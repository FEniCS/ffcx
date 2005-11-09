__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2004-11-03"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

# FFC common modules
from ffc.common.debug import *
from ffc.common.util import *

# FFC compiler modules
from algebra import *
from reassign import *
from multiindex import *

class GeometryTensor:

    """A GeometryTensor represents the geometry tensor of a
    multi-linear form computed on the current cell.

        a            - secondary multiindex
        b            - auxiliary multiindex
        rank         - rank of the tensor
        constants    - a list of Constants
        coefficients - a list of Coefficients
        transforms   - a list of Transforms"""

    def __init__(self, product):
        "Create GeometryTensor."

        # Check that we get a Product
        if not isinstance(product, Product):
            raise RuntimeError, "GeometryTensor must be created from Product."

        # Get data from Product
        self.constants = listcopy(product.constants)
        self.coefficients = listcopy(product.coefficients)
        self.transforms = listcopy(product.transforms)

        # Create MultiIndices
        self.a = self.__create_index("secondary")
        self.b = self.__create_index("geometry tensor auxiliary")

        # Get rank
        self.rank = self.a.rank

        debug("Created geometry tensor: a%s, b%s" % \
              (str(self.a.dims), str(self.b.dims)), 1)

        return

    def __create_index(self, type):
        "Find dimensions and create MultiIndex."
        # Compute rank
        rank = max([max_index(c, type) for c in self.coefficients] + \
                   [max_index(t, type) for t in self.transforms] + [-1]) + 1
        
        # Compute all dimensions
        dims = [self.__find_dim(i, type) for i in range(rank)]
        # Create MultiIndex
        return MultiIndex(dims)

    def __find_dim(self, i, type):
        "Find dimension of given Index."
        index = Index(0)
        index.index = i
        index.type = type
        # Check coefficients
        for c in self.coefficients:
            if c.index == index:
                return c.e1.spacedim()
        # Check transforms
        for t in self.transforms:
            if t.index0 == index or t.index1 == index:
                return t.element.shapedim()
        # Didn't find dimension
        raise RuntimeError, "Unable to find dimension for Index " + str(index)

    def __call__(self, a, format):
        "Return given element of geometry tensor."
        aux = "geometry tensor auxiliary"
        # Compute product of factors outside sum
        factors = []
        for c in self.constants:
            if c.inverted:
                factors += ["(1.0/" + format.format["constant"](c.number.index) + ")"]
            else:
                factors += [format.format["constant"](c.number.index)]
        for c in self.coefficients:
            if not c.index.type == aux:
                factors += [format.format["coefficient"](c.n1.index, c.index([], a, [], []))]
        for t in self.transforms:
            if not (t.index0.type == aux or  t.index1.type == aux):
                factors += [format.format["transform"](t.index0([], a, [], []), \
                                                       t.index1([], a, [], []))]
        product = format.format["multiplication"](factors)
        if product: f0 = [product]
        else: f0 = []
        # Compute sum of products inside sum
        terms = []
        for b in self.b.indices:
            factors = []
            for c in self.coefficients:
                if c.index.type == aux:
                    factors += [format.format["coefficient"](c.n1.index, c.index([], a, [], b))]
            for t in self.transforms:
                if t.index0.type == aux or t.index1.type == aux:
                    factors += [format.format["transform"](t.index0([], a, [], b), \
                                                           t.index1([], a, [], b))]
            terms += [format.format["multiplication"](factors)]
        sum = format.format["sum"](terms)
        if sum: sum = format.format["grouping"](sum)
        if sum: f1 = [sum]
        else: f1 = []
        # Compute product of all factors
        return format.format["multiplication"]([f for f in [format.format["determinant"]] + f0 + f1])

    def __repr__(self):
        "Print nicely formatted representation of GeometryTensor."
        c = "".join([c.__repr__() for c in self.coefficients])
        t = "".join([t.__repr__() for t in self.transforms])
        return c + t
