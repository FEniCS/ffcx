__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2004-11-04"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

class Integrand:
    
    """This class wraps around a list of BasisFunctions to create a
    callable object."""

    def __init__(self, basisfunctions, indices0, indices1, vscaling, dscaling):
        "Create Integrand."
        self.basisfunctions = basisfunctions
        self.indices0 = indices0
        self.indices1 = indices1
        self.vscaling = vscaling
        self.dscaling = dscaling
        return

    def __call__(self, x):
        "Evaluate product at given point."
        tmp = 1.0
        for v in self.basisfunctions:
            tmp = tmp * self.__eval(v, x)
        return tmp * self.vscaling

    def __eval(self, v, x):
        "Evaluate BasisFunction at given point."

        # Reset derivative scaling
        scaling = 1.0

        # Get basis index of BasisFunction
        index = v.index(self.indices0, self.indices1, [], [])

        # Get component index of BasisFunction
        component = [i(self.indices0, self.indices1, [], []) for i in v.component]

        # Get basis function
        if len(component) > 0:
            w = v.element.basis[index][component]
        else:
            w = v.element.basis[index]

        # Differentiate the basis function
        for d in v.derivatives:
            i = d.index(self.indices0, self.indices1, [], [])
            w = w.deriv(i)
            scaling *= self.dscaling

        # Evaluate basis function
        return w(x) * scaling
