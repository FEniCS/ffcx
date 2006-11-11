__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2004-11-03 -- 2006-05-07"
__copyright__ = "Copyright (C) 2004-2006 Anders Logg"
__license__  = "GNU GPL Version 2"

# Modified by Garth N. Wells 2006

# Python modules
import numpy
import time

# FFC common modules
from ffc.common.exceptions import *
from ffc.common.progress import *
from ffc.common.debug import *
from ffc.common.util import *

# FFC compiler modules
from algebra import *
from reassign import *
from multiindex import *
from integrator import *
from monomialintegration import *

class ReferenceTensor:

    """A ReferenceTensor represents the reference tensor of a
    multi-linear form computed on the reference cell.

    A ReferenceTensor holds the following data:

        i              - primary multiindex
        a              - secondary multiindex
        b              - auxiliary multiindex
        rank           - rank of the tensor
        A0             - the precomputed reference tensor
        numeric        - a numeric constant (float)
        basisfunctions - a list of BasisFunctions
        integral       - an Integral
        cputime        - time to compute the reference tensor"""

    def __init__(self, product, facet):
        "Create ReferenceTensor."

        # Check that we get a Product
        if not isinstance(product, Product):
            raise RuntimError, "ReferenceTensor must be created from Product."

        # Check that the Product contains an Integral
        if product.integral == None:
            raise FormError, (product, "Missing integral in term.")

        # Get data from Product
        self.numeric = product.numeric
        self.basisfunctions = listcopy(product.basisfunctions)
        self.integral = product.integral

        # Create MultiIndices
        self.i = self.__create_index("primary")
        self.a = self.__create_index("secondary")
        self.b = self.__create_index("reference tensor auxiliary")

        # Compute reference tensor
        t = time.time()
        self.A0 = integrate(product, facet)

        # Report time to compute the reference tensor
        self.cputime = time.time() - t
        debug("Reference tensor computed in %.3g seconds" % self.cputime)

        # For testing
        #B = self.__compute_reference_tensor()
        #print "Maximum error: %.3e" % max(abs(numpy.ravel(self.A0 - B)))

        # Compute reference tensor (old version)
        #self.A0 = self.__compute_reference_tensor()

        # Get rank
        self.rank = self.i.rank + self.a.rank

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
                return v.element.spacedim()
            # Check component Indices
            for j in range(len(v.component)):
                if v.component[j] == index:
                    return v.element.tensordim(j)
            # Check Derivatives
            for d in v.derivatives:
                if d.index == index:
                    return d.element.shapedim()
        # Didn't find dimension
        raise RuntimeError, "Unable to find dimension for Index " + str(index)

    def __compute_reference_tensor(self):
        "Compute reference tensor."

        # This is the old (slow) version, not used but still here for testing

        # Make sure that the iteration is not empty
        iindices = self.i.indices or [[]]
        aindices = self.a.indices or [[]]
        bindices = self.b.indices or [[]]

        # Create tensor
        A0 = numpy.zeros(self.i.dims + self.a.dims, dtype = numpy.float)
        debug("Number of entries in reference tensor: %d" % numpy.size(A0), 1)

        # Create quadrature rule
        integrate = Integrator(self.basisfunctions)

        # Count the number of integrals
        n = numpy.product(self.i.dims) * \
            numpy.product(self.a.dims) * \
            numpy.product(self.b.dims)
        debug("Computing %d integrals, this may take some time" % n)

        # Iterate over all combinations of indices
        debug("Computing reference tensor", 2)
        progress = Progress(n)
        for i in iindices:
            debug("i = " + str(i), 3)
            for a in aindices:
                debug("  a = " + str(a), 3)
                integral = 0.0
                for b in bindices:
                    debug("    b = " + str(b), 3)
                    integral += integrate(self.basisfunctions, i, a, b)
                    progress += 1
                integral *= self.numeric
                A0[i + a] = integral
                debug("  integral = " + str(integral), 3)
        return A0

    def  __call__(self, i, a = []):
        "Return given element of reference tensor."
        return self.A0[i + a]

    def __repr__(self):
        "Print nicely formatted representation of ReferenceTensor."
        v = "*".join([v.__repr__() for v in self.basisfunctions])
        return v + "*" + self.integral.__repr__()
