__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2004-11-04"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

# FFC modules
from tensorspace import ZeroFunction

class Integrand:
    
    """This class wraps around a list of BasisFunctions to create a
    callable object."""

    def __init__(self, basisfunctions, iindices, aindices, bindices, vscaling, dscaling):
        "Create Integrand."
        self.basisfunctions = basisfunctions
        self.iindices = iindices
        self.aindices = aindices
        self.bindices = bindices
        self.vscaling = vscaling
        self.dscaling = dscaling
        return

    def __call__(self, x):
        "Evaluate product at given point."
        tmp = 1.0
        for v in self.basisfunctions:
            tmp = tmp * self.__eval(v, x)
        return tmp * self.vscaling

    def iszero(self):
        "Return true iff integrand is zero."
        for v in self.basisfunctions:
            if isinstance(self.__getcomponent(v), ZeroFunction):
                return True
        return False

    def __eval(self, v, x):
        "Evaluate BasisFunction at given point."

        # Reset derivative scaling
        scaling = 1.0

        # Get current component of current basis function
        w = self.__getcomponent(v)

        # Differentiate the basis function
        for d in v.derivatives:
            i = d.index(self.iindices, self.aindices, self.bindices, [])
            w = w.deriv(i)
            scaling *= self.dscaling

        # Evaluate basis function
        return w(x) * scaling
    
    def __getcomponent(self, v):
        "Return current component of current basis function."
        
        # Get basis index of BasisFunction
        index = v.index(self.iindices, self.aindices, self.bindices, [])

        # Get component index of BasisFunction
        component = [i(self.iindices, self.aindices, self.bindices, []) for i in v.component]
        
        # Get component of BasisFunction
        if len(component) > 0:
            # Tensor case
            w = v.element.basis[index]
            for i in component:
                w = w[i]
        else:
            # Scalar case
            w = v.element.basis[index]

        return w
