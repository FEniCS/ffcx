"An algebra for multi-linear forms."

__author__ = "Anders Logg (logg@tti-c.org)"
__version__ = "0.0.1"
__date__ = "2004-09-27"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

from derivative import Derivative
from index import Index

# FFC common modules
from ffc.common.log import error

class Element:
    "Base class for elements of the algebra."
    pass

class BasisFunction(Element):
    """A BasisFunction represents a basis function on the reference
    cell and can be either an argument in the multi-linear form or an
    auxiliary basis function that will be removed in the
    contraction."""

    def __init__(self, index = None):
        if index == None:
            self.index = Index("primary")
        else:
            self.index = Index(index)
        return

    def __add__(self, other):
        "Operator: BasisFunction + Element"
        if isinstance(other, BasisFunction):
            w = Sum()
            w.products = [Product(self), Product(other)]
            return w
        elif isinstance(other, Factor):
            w = Sum()
            w.products = [Product(self), Product(other)]
            return w
        elif isinstance(other, Product):
            w = Sum()
            w.products = [Product(self), other]
            return w
        elif isinstance(other, Sum):
            w = Sum()
            w.products = [Product(self)] + other.products
            return w
        else:
            error("Can't add BasisFunction with " + str(other))
        return

    def __sub__(self, other):
        "Operator: BasisFunction - Element"
        # FIXME: Remember to modify geometry tensor with a -
        if isinstance(other, BasisFunction):
            w = Sum()
            w.products = [Product(self), Product(other)]
            return w
        elif isinstance(other, Factor):
            w = Sum()
            w.products = [Product(self), Product(other)]
            return w
        elif isinstance(other, Product):
            w = Sum()
            w.products = [Product(self), other]
            return w
        elif isinstance(other, Sum):
            w = Sum()
            w.products = [Product(self)] + other.products
            return w
        else:
            error("Can't subtract BasisFunction with " + str(other))
        return

    def __mul__(self, other):
        "Operator: BasisFunction * Element"
        if isinstance(other, BasisFunction):
            w = Product()
            w.factors = [Factor(self), Factor(other)]
            return w
        elif isinstance(other, Factor):
            w = Product()
            w.factors = [Factor(self), other]
            return w
        elif isinstance(other, Product):
            w = Product()
            w.factors = [Factor(self)] + other.factors
            return w
        elif isinstance(other, Sum):
            w = Sum()
            for p in other.products:
                w.products += [self*p]
            return w
        else:
            error("Can't multiply BasisFunction with " + str(other))
        return

    def dx(self, index = None):
        "Operator: (d/dx)BasisFunction in given coordinate direction."
        w = Factor(self)
        w.derivatives = [Derivative(index)]
        return w

    def __repr__(self):
        "Print nicely formatted representation of BasisFunction."
        if self.index.type == "primary":
            return "v" + str(self.index)
        else:
            return "w" + str(self.index)
        
class Factor(Element):
    """A Factor represents a (possibly) differentiated BasisFunction
    on the reference cell."""
    
    def __init__(self, other = None):
        if other == None:
            self.basisfunction = None
            self.derivatives = []
        elif isinstance(other, BasisFunction):
            self.basisfunction = other
            self.derivatives = []
        elif  isinstance(other, Factor):
            self.basisfunction = other.basisfunction
            self.derivatives = other.derivatives
        else:
            error("Unable to create Factor from " + str(other))
        return

    def __add__(self, other):
        if isinstance(other, BasisFunction):
            w = Sum()
            w.products = [Product(self), Product(other)]
            return w
        elif isinstance(other, Factor):
            w = Sum()
            w.products = [Product(self), Product(other)]
            return w
        elif isinstance(other, Product):
            w = Sum()
            w.products = [Product(self), other]
            return w
        elif isinstance(other, Sum):
            w = Sum()
            w.products = [Product(self)] + other.products
            return w
        else:
            error("Can't add Factor with " + str(other))
        return

    def __sub__(self, other):
        # FIXME: Remember to modify geometry tensor with a -
        if isinstance(other, BasisFunction):
            w = Sum()
            w.products = [Product(self), Product(other)]
            return w
        elif isinstance(other, Factor):
            w = Sum()
            w.products = [Product(self), Product(other)]
            return w
        elif isinstance(other, Product):
            w = Sum()
            w.products = [Product(self), other]
            return w
        elif isinstance(other, Sum):
            w = Sum()
            w.products = [Product(self)] + other.products
            return w
        else:
            error("Can't add Factor with " + str(other))
        return
    
    def __mul__(self, other):
        if isinstance(other, BasisFunction):
            w = Product()
            w.factors = [self, Factor(other)]
            return w
        elif isinstance(other, Factor):
            w = Product()
            w.factors = [self, other]
            return w
        elif isinstance(other, Product):
            w = Product()
            w.factors = [self] + other.factors
            return w
        elif isinstance(other, Sum):
            w = Sum()
            for p in other.products:
                w.products += [self*p]
            return w
        else:
            error("Can't multiply Factor with " + str(other))
        return

    def dx(self, index = None):
        "Operator: (d/dx)Factor in given coordinate direction."
        w = Factor(self)
        w.derivatives = [Derivative(index)] + w.derivatives
        return w
        
    def __repr__(self):
        "Print nicely formatted representation of Factor."
        output = ""
        if len(self.derivatives) > 0:
            output += "("
            for d in self.derivatives:
                output += d.__repr__()
            output += self.basisfunction.__repr__()
            output += ")"
        else:
            output += self.basisfunction.__repr__()
        return output
    
class Product(Element):
    "A Product represents a product of Factors."

    def __init__(self, other = None):
        if other == None:
            self.factors = []
        elif isinstance(other, BasisFunction):
            self.factors = [Factor(other)]
        elif isinstance(other, Factor):
            self.factors = [other]
        elif isinstance(other, Product):
            self.factors = other.factors
        else:
            error("Unable to create Product from " + str(other))
        return
    
    def __add__(self, other):
        if isinstance(other, BasisFunction):
            w = Sum()
            w.products = [self, Product(other)]
            return w
        elif isinstance(other, Factor):
            w = Sum()
            w.products = [self, Product(other)]
            return w
        elif isinstance(other, Product):
            w = Sum()
            w.products = [self, other]
            return w
        elif isinstance(other, Sum):
            w = Sum()
            w.products = [self] + other.products
            return w
        else:
            error("Can't add Product with " + str(other))
        return

    def __sub__(self, other):
        # FIXME: Remember to modify geometry tensor with a -
        if isinstance(other, BasisFunction):
            w = Sum()
            w.products = [self, Product(other)]
            return w
        elif isinstance(other, Factor):
            w = Sum()
            w.products = [self, Product(other)]
            return w
        elif isinstance(other, Product):
            w = Sum()
            w.products = [self, other]
            return w
        elif isinstance(other, Sum):
            w = Sum()
            w.products = [self] + other.products
            return w
        else:
            error("Can't subtract Product with " + str(other))
        return
    
    def __mul__(self, other):
        if isinstance(other, BasisFunction):
            w = Product()
            w.factors = self.factors + [Factor(other)]
            return w
        elif isinstance(other, Factor):
            w = Product()
            w.factors = self.factors + [other]
            return w
        elif isinstance(other, Product):
            w = Product()
            w.factors = self.factors + other.factors
            return w
        elif isinstance(other, Sum):
            w = Sum()
            for p in other.products:
                w.products += [self*p]
            return w
        else:
            error("Can't multiply Product with " + str(other))
        return

    def dx(self, index = None):
        "Operator: (d/dx)Product in given coordinate direction."
        w = Sum()
        for i in range(len(self.factors)):
            p = Product()
            for j in range(len(self.factors)):
                if i == j:
                    p.factors += [self.factors[i].dx(index)]
                else:
                    p.factors += [self.factors[j]]
            w.products += [p]
        return w

    def rank(self):
        "Return rank [r0, r1, rtot] of tensor represented by the Product."
        
        

    def __repr__(self):
        "Print nicely formatted representation of Product."
        output = ""
        for i in range(len(self.factors)):
            if i < (len(self.factors) - 1):
                output += (self.factors[i].__repr__() + "*")
            else:
                output += self.factors[i].__repr__()
        return output

class Sum(Element):
    "A Sum represents a sum of Products."
    
    def __init__(self, other = None):
        if other == None:
            self.products = []
        elif isinstance(other, BasisFunction):
            self.products = [Product(other)]
        elif isinstance(other, Factor):
            self.products = [Product(other)]
        elif isinstance(other, Product):
            self.products = [other]
        elif isinstance(other, Sum):
            self.products = other.products
        else:
            error("Unable to create Sum from " + str(other))
        return

    def __add__(self, other):
        if isinstance(other, BasisFunction):
            w = Sum()
            w.products = self.products + [Product(other)]
            return w
        elif isinstance(other, Factor):
            w = Sum()
            w.products = self.products + [Product(other)]
            return w
        elif isinstance(other, Product):
            w = Sum()
            w.products = self.products + [other]
            return w
        elif isinstance(other, Sum):
            w = Sum()
            w.products = self.products + other.products
            return w
        else:
            error("Can't add Sum with " + str(other))
        return

    def __sub__(self, other):
        # FIXME: Remember to modify geometry tensor with a -
        if isinstance(other, BasisFunction):
            w = Sum()
            w.products = self.products + [Product(other)]
            return w
        elif isinstance(other, Factor):
            w = Sum()
            w.products = self.products + [Product(other)]
            return w
        elif isinstance(other, Product):
            w = Sum()
            w.products = self.products + [other]
            return w
        elif isinstance(other, Sum):
            w = Sum()
            w.products = self.products + other.products
            return w
        else:
            error("Can't add Sum with " + str(other))
        return
    
    def __mul__(self, other):
        if isinstance(other, BasisFunction):
            w = Sum()
            for p in self.products:
                w.products += [p*other]
            return w
        elif isinstance(other, Factor):
            w = Sum()
            for p in self.products:
                w.products += [p*other]
            return w
        elif isinstance(other, Product):
            w = Sum()
            for p in self.products:
                w.products += [p*other]
            return w
        elif isinstance(other, Sum):
            w = Sum()
            for p in self.products:
                for q in other.products:
                    w.products += [p*q]
            return w
        else:
            error("Can't multiply Sum with " + str(other))
        return

    def dx(self, index = None):
        "Operator: (d/dx)Sum in given coordinate direction."
        w = Sum()
        for p in self.products:
            w += p.dx(index)
        return w

    def __repr__(self):
        "Print nicely formatted representation of Sum."
        output = ""
        for i in range(len(self.products)):
            if i < (len(self.products) - 1):
                output += (self.products[i].__repr__() + " + ")
            else:
                output += self.products[i].__repr__()
        return output

if __name__ == "__main__":

    print "Testing algebra"
    print "---------------"

    u = BasisFunction()
    v = BasisFunction()

    print "Testing long expression:"
    w = (((u*(u*v + u)*u + u*v)*u.dx(0)).dx(2)).dx(1)
    print w

    print

    print "Testing derivative of product:"
    w = (u*v*u).dx(0)
    print w

    print

    print "Testing derivative of sum:"
    w = (u + v).dx(0)
    print w

    print
    
    print "Testing Poisson:"
    i = Index()
    w = u.dx(i)*v.dx(i)
    print w
