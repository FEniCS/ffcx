__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2004-11-03"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

# Python modules
from Numeric import *

# FFC common modules
from ffc.common.progress import *
from ffc.common.debug import *
from ffc.common.util import *

# FFC compiler modules
from algebra import *
from integral import *
from reassign import *
from multiindex import *
from integrator import *

class ReferenceTensor:

    """A ReferenceTensor represents the reference tensor of a
    multi-linear form computed on the reference cell.

    A ReferenceTensor holds the following data:

        i              - primary multiindex
        a              - secondary multiindex
        b              - auxiliary multiindex
        A0             - the precomputed reference tensor
        numeric        - a numeric constant (float)
        basisfunctions - a list of BasisFunctions
        integral       - an Integral"""

    def __init__(self, product):
        "Create ReferenceTensor."

        # Check that we get a Product
        if not isinstance(product, Product):
            raise RuntimeError, "ReferenceTensor must be created from Product."

        # Check that the Product contains an Integral
        if product.integral == None:
            raise RuntimeError, "Missing integral."

        # Get data from Product
        self.numeric = product.numeric
        self.basisfunctions = listcopy(product.basisfunctions)
        self.integral = product.integral

        # Create MultiIndices
        self.i = self.__create_index("primary")
        self.a = self.__create_index("secondary")
        self.b = self.__create_index("reference tensor auxiliary")

        # Compute ReferenceTensor
        self.A0 = self.__compute_reference_tensor()

        debug("Created reference tensor: i%s, a%s, b%s" % \
              (str(self.i.dims), str(self.a.dims), str(self.b.dims)), 1)
        
        return

    def __create_index(self, type):
        "Find dimensions and create MultiIndex."
        # Compute rank
        rank = max([max_index(v, type) for v in self.basisfunctions] + [-1]) + 1
        # Compute all dimensions
        dims = [self.__find_dim(i, type) for i in range(rank)]
        # Create MultiIndex
        return MultiIndex(dims)

    def __find_dim(self, i, type):
        "Find dimension of given Index."
        index = Index(i)
        index.type = type
        for v in self.basisfunctions:
            # Check basis Index
            if v.index == index:
                return v.element.spacedim
            # Check component Indices
            for j in range(len(v.component)):
                if v.component[j] == index:
                    return v.element.tensordims[j]
            # Check Derivatives
            for d in v.derivatives:
                if d.index == index:
                    return d.element.shapedim
        # Didn't find dimension
        raise RuntimeError, "Unable to find dimension for Index " + str(index)

    def __compute_reference_tensor(self):
        "Compute reference tensor."

        # Make sure that the iteration is not empty
        iindices = self.i.indices
        aindices = self.a.indices or [[]]
        bindices = self.b.indices or [[]]

        # Create tensor
        A0 = zeros(self.i.dims + self.a.dims, Float)

        # Create quadrature rule
        integrate = Integrator(self.basisfunctions)

        # Count the number of integrals
        n = product(self.i.dims) * product(self.a.dims) * product(self.b.dims)
        debug("Computing %d integrals, this may take some time" % n)

        # Iterate over all combinations of indices
        debug("Computing reference tensor", 2)
        progress = Progress(n)
        for i in iindices:
            debug("i = " + str(i), 2)
            for a in aindices:
                debug("  a = " + str(a), 2)
                integral = 0.0
                for b in bindices:
                    debug("    b = " + str(b), 2)
                    integral += integrate(self.basisfunctions, i, a, b)
                    progress += 1
                integral *= self.numeric
                A0[i + a] = integral
                debug("  integral = " + str(integral), 2)
        return A0

    def  __call__(self, i, a = []):
        "Return given element of reference tensor."
        return self.A0[i + a]

    def __repr__(self):
        "Print nicely formatted representation of ReferenceTensor."
        v = "*".join([v.__repr__() for v in self.basisfunctions])
        return v + "*" + self.integral.__repr__()
