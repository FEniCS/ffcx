__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2004-11-16"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

# FIAT modules
from FIAT import shapes

class ZeroFunction:

    "A function that is always zero."

    def __call__(self, x):
        "Evaluate function."
        return 0.0

    def deriv(self, i):
        "Differentiate function."
        return self

class TensorSpace:

    """This is a temporary implementation of a general tensor-valued
    finite element space until general tensor-valued elements are
    added to FIAT."""

    def __init__(self, scalarbasis, dims):
        "Create TensorSpace."

        # Only vector-valued elements implemented so far
        if not len(dims) == 1:
            raise RuntimeError, "General tensor-valued elements not implemented."
        dim = dims[0]

        # Create vector-valued basis
        self.basis = []
        for v in scalarbasis:
            for i in range(dim):
                w = [ZeroFunction() for i in range(dim)]
                w[i] = v
                self.basis += [w]

        return

    def getBasis(self):
        return self.basis
