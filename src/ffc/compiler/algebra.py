"""An algebra for multi-linear forms. Objects of the following classes
are elements of the algebra:

    BasisFunction - basic building block
    Function      - linear combination of BasisFunctions
    Product       - product of BasisFunctions
    Sum           - sum of Products

Each element of the algebra, BasisFunction, Function, Product, and
Sum, can be either scalar or tensor-valued. The following operations
are supported for all elements of the algebra:

    Binary +      (tensor ranks of operands must match)
    Binary -      (tensor ranks of operands must match)
    Binary *      (both operands must be scalar)
    Unary  +      (operand scalar or tensor-valued)
    Unary  -      (operand scalar or tensor-valued)
    Unary  []     (operand must be tensor-valued)
    Unary  d/dx   (operand scalar or tensor-valued)"""

__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2004-09-27"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

# FFC modules
from finiteelement import FiniteElement
from coefficient import Coefficient
from derivative import Derivative
from transform import Transform
from integral import Integral
from index import Index

class Element:

    "Base class for elements of the algebra."

    def __add__(self, other):
        "Operator: Element + Element"
        return Sum(self) + other

    def __radd__(self, other):
        "Operator: Element + Element (other + self)"
        return Sum(self) + other

    def __sub__(self, other):
        "Operator: Element - Element"
        return Sum(self) + (-other)

    def __rsub__(self, other):
        "Operator: Element + Element (other - self)"
        return Sum(-self) + other

    def __mul__(self, other):
        "Operator: Element * Element"
        return Sum(self) * other

    def __rmul__(self, other):
        "Operator: Element * Element (other * self)"
        return Sum(self) * other

    def __div__(self, other):
        "Operator: Element / Element (only works if other is scalar)."
        return Sum(self) * (1.0/other)

    def __pos__(self):
        "Operator: +Element"
        return Sum(self)

    def __neg__(self):
        "Operator: -Element"
        return -Sum(self)

    def  __getitem__(self, component):
        "Operator: Element[component], pick given component."
        return Sum(self)[component]

    def dx(self, index = None):
        "Operator: (d/dx)Element in given coordinate direction."
        return Sum(self).dx(index)

    def __repr__(self):
        "Print nicely formatted representation of Element."
        return Sum(self).__repr__()

    def rank(self):
        "Return value rank of Element."
        return Sum(self).rank()

class Function(Element):

    """A Function represents a projection of a given function onto a
    finite element space, expressed as a linear combination of
    BasisFunctions.

    A Function holds the following data:

        element - a FiniteElement
        name    - a string identifying the function."""

    def __init__(self, element, name = None):
        "Create Function."
        if isinstance(element, Function):
            # Create Function from Function (copy constructor)
            self.element = element.element
            self.name = element.name
        else:
            # Create Function for given FiniteElement with given name
            self.element = element
            self.name = element.name
        return

class BasisFunction(Element):

    """A BasisFunction represents a possibly differentiated component
    of a basis function on the reference cell.

    A BasisFunction holds the following data:

        element     - a FiniteElement
        index       - a basis Index
        component   - a list of component Indices
        derivatives - a list of Derivatives"""

    def __init__(self, element, index = None):
        "Create BasisFunction."
        if index == None and isinstance(element, BasisFunction):
            # Create BasisFunction from BasisFunction (copy constructor)
            self.element = element.element
            self.index = Index(element.index)
            self.component = [] + element.component
            self.derivatives = [] + element.derivatives
        elif index == None:
            # Create BasisFunction with primary Index (default)
            self.element = element
            self.index = Index("primary")
            self.component = []
            self.derivatives = []
        else:
            # Create BasisFunction with specified Index
            self.element = element
            self.index = Index(index)
            self.component = []
            self.derivatives = []
        return

    def  __getitem__(self, component):
        "Operator: BasisFunction[component], pick given component."
        rank = self.element.rank
        if self.component or rank == 0:
            raise RuntimeError, "Cannot pick component of scalar BasisFunction."
        w = BasisFunction(self)
        if isinstance(component, list):
            if not rank == len(component):
                raise RuntimeError, "Illegal component index, does not match rank."
            w.component = [] + component
        else:
            if not rank == 1:
                raise RuntimeError, "Illegal component index, does not match rank."
            w.component = [Index(component)]
        return w

    def dx(self, index = None):
        "Operator: (d/dx)BasisFunction in given coordinate direction."
        i = Index() # Create new secondary index
        w = Product(self)
        w.basisfunctions[0].derivatives.insert(0, Derivative(self.element, i))
        w.transforms.insert(0, Transform(self.element, i, index))
        return w

    def __repr__(self):
        "Print nicely formatted representation of BasisFunction."
        d = "".join([d.__repr__() for d in self.derivatives])
        i = self.index.__repr__()
        if self.component:
            c = "[" + ",".join([c.__repr__() for c in self.component]) + "]"
        else:
            c = ""
        if len(self.derivatives) > 0:
            return "(" + d + "v" + i + c + ")"
        else:
            return d + "v" + i + c

    def rank(self):
        "Return value rank of BasisFunction."
        if self.component:
            return 0
        else:
            return self.element.rank

    def indexcall(self, foo, args = None):
        "Call given function on all Indices."
        self.index.indexcall(foo, args)
        [i.indexcall(foo, args) for i in self.component]
        [d.indexcall(foo, args) for d in self.derivatives]
        return

class Product(Element):

    """A Product represents a product of factors, including
    BasisFunctions and Functions.

    A Product holds the following data:

        constant       - a numeric constant (float)
        coefficients   - a list of Coefficients
        transforms     - a list of Transforms
        basisfunctions - a list of BasisFunctions
        integral       - an Integral"""

    def __init__(self, other = None):
        "Create Product."
        if other == None:
            # Create default Product (unity)
            self.constant = 1.0
            self.coefficients = []
            self.transforms = []
            self.basisfunctions = []
            self.integral = None
        elif isinstance(other, int) or isinstance(other, float):
            # Create Product from scalar
            self.constant = float(other)
            self.coefficients = []
            self.transforms = []
            self.basisfunctions = []
            self.integral = None
        elif isinstance(other, Function):
            # Create Product from Function
            index = Index()
            self.constant = 1.0
            self.coefficients = [Coefficient(other.element, other.name, index)]
            self.transforms = []
            self.basisfunctions = [BasisFunction(other.element, index)]
            self.integral = None
        elif isinstance(other, BasisFunction):
            # Create Product from BasisFunction
            self.constant = 1.0
            self.coefficients = []
            self.transforms = []
            self.basisfunctions = [BasisFunction(other)]
            self.integral = None
        elif isinstance(other, Product):
            # Create Product from Product (copy constructor)
            self.constant = other.constant
            self.coefficients = [] + other.coefficients
            self.transforms = [] + other.transforms
            self.basisfunctions = [] + other.basisfunctions
            self.integral = other.integral
        else:
            raise RuntimeError, "Unable to create Product from " + str(other)
        return
    
    def __mul__(self, other):
        "Operator: Product * Element"
        if isinstance(other, Integral):
            if not self.integral == None:
                raise RuntimeError, "Integrand can only be integrated once."
            w = Product(self)
            w.integral = other
            return w
        elif isinstance(other, Sum):
            return Sum(self) * Sum(other)
        else:
            w0 = Product(self)
            w1 = Product(other)
            if not w0.rank() == w1.rank() == 0:
                raise RuntimeError, "Illegal ranks for product, must be scalar."
            if w0.integral and w1.integral:
                raise RuntimeError, "Integrand can only be integrated once."
            w = Product()
            w.constant = w0.constant * w1.constant
            w.coefficients = w0.coefficients + w1.coefficients
            w.transforms = w0.transforms + w1.transforms
            w.basisfunctions = w0.basisfunctions + w1.basisfunctions
            w.integral = w0.integral or w1.integral
            return w

    def __neg__(self):
        "Operator: -Product"
        w = Product(self)
        w.constant = -w.constant
        return w

    def  __getitem__(self, component):
        "Operator: Product[component], pick given component."
        if not len(self.basisfunctions) == 1:
            # Always scalar if product of more than one basis function
            raise RuntimeError, "Cannot pick component of scalar Product."
        w = Product(self)
        w.basisfunctions[0] = w.basisfunctions[0][component]
        return w

    def dx(self, index = None):
        "Operator: (d/dx)Product in given coordinate direction."
        w = Sum()
        for i in range(len(self.basisfunctions)):
            p = Product(self)
            p.basisfunctions = []
            for j in range(len(self.basisfunctions)):
                if i == j:
                    p = p * self.basisfunctions[i].dx(index)
                else:
                    p = p * self.basisfunctions[j]
            w = w + p
        return w

    def __repr__(self):
        "Print nicely formatted representation of Product."
        if not (self.coefficients or self.transforms or self.basisfunctions):
            return str(self.constant)
        if self.constant == -1.0:
            s = "-"
        elif not self.constant == 1.0:
            s = str(self.constant)
        else:
            s = ""
        c = "".join([c.__repr__() for c in self.coefficients])
        t = "".join([t.__repr__() for t in self.transforms])
        v = "*".join([v.__repr__() for v in self.basisfunctions])
        if not self.integral == None:
            i = "*" + self.integral.__repr__()
        else:
            i = ""
        return s + c + t + v + i

    def rank(self):
        "Return value rank of Product."
        if not self.basisfunctions:
            return 0
        if len(self.basisfunctions) > 1:
            for v in self.basisfunctions:
                if not v.rank() == 0:
                    raise RuntimeError, "Illegal rank for BasisFunction of Product (non-scalar)."
        return self.basisfunctions[0].rank()
            
    def indexcall(self, foo, args = None):
        "Call given function on all Indices."
        [c.indexcall(foo, args) for c in self.coefficients]
        [t.indexcall(foo, args) for t in self.transforms]
        [v.indexcall(foo, args) for v in self.basisfunctions]
        return

class Sum(Element):

    """A Sum represents a sum of Products. Each Product will be
    compiled separately, since different Products are probably of
    different rank.

    A Sum holds the following data:

        products - a list of Products (terms) in the Sum"""
    
    def __init__(self, other = None):
        "Create Sum."
        if other == None:
            # Create default Sum (zero)
            self.products = []
        elif isinstance(other, int) or isinstance(other, float):
            # Create Sum from scalar
            self.products = [Product(other)]
        elif isinstance(other, BasisFunction):
            # Create Sum from BasisFunction
            self.products = [Product(other)]
        elif isinstance(other, Function):
            # Create Sum from Function
            self.products = [Product(other)]
        elif isinstance(other, Product):
            # Create Sum from Product
            self.products = [Product(other)]
        elif isinstance(other, Sum):
            # Create Sum from Sum (copy constructor)
            self.products = [] + other.products
        else:
            raise RuntimeError, "Unable to create Sum from " + str(other)
        return

    def __add__(self, other):
        "Operator: Sum + Element"
        w0 = Sum(self)
        w1 = Sum(other)
        if not w0.rank() == w1.rank():
            raise RuntimeError, "Non-matching ranks."
        w = Sum()
        w.products = w0.products + w1.products
        return w
    
    def __mul__(self, other):
        "Operator: Sum * Element"
        if isinstance(other, Sum):
            w = Sum()
            w.products = [p*q for p in self.products for q in other.products]
            return w
        else:
            w = Sum()
            w.products = [p*other for p in self.products]
            return w

    def __neg__(self):
        "Operator: -Sum"
        w = Sum(self)
        w.products = [-p for p in w.products]
        return w

    def  __getitem__(self, component):
        "Operator: indexing, pick given component."
        w = Sum(self)
        w.products = [p[component] for p in w.products]
        return w

    def dx(self, index = None):
        "Operator: (d/dx)Sum in given coordinate direction."
        w = Sum()
        w.products = [p.dx(index) for p in self.products]
        return w

    def __repr__(self):
        "Print nicely formatted representation of Sum."
        return " + ".join([p.__repr__() for p in self.products])

    def rank(self):
        "Return value rank of Sum."
        if not self.products:
            return 0
        for j in range(len(self.products) - 1):
            if not self.products[j].rank() == self.products[j + 1].rank():
                raise RuntimeError, "Non-matching ranks."
        return self.products[0].rank()
  
    def indexcall(self, foo, args = None):
        "Call given function on all Indices."
        [p.indexcall(foo, args) for p in self.products]
        return

if __name__ == "__main__":

    print "Testing algebra"
    print "---------------"

    element = FiniteElement("Lagrange", "triangle", 1)
    u = BasisFunction(element)
    v = BasisFunction(element)
    dx = Integral("interior")
    ds = Integral("boundary")

    print "Testing long expression:"
    w = 2*(u*(u/2 + v)*u + u*v)*u.dx(0).dx(2).dx(1)
    print w

    print

    print "Testing derivative of product:"
    w = (u*v*u).dx(0)*3
    print w

    print

    print "Testing derivative of sum:"
    w = (u + v).dx(0)
    print w

    print
    
    print "Testing Poisson:"
    i = Index()
    w = u.dx(i)*v.dx(i)*dx + u*v*ds
    print w

    print

    print "Testing Biharmonic:"
    i = Index()
    j = Index()
    w = u.dx(i).dx(i)*v.dx(j).dx(j)
    print w

    print

    print "Testing Poisson system:"
    element = FiniteElement("Lagrange", "triangle", 1, 3)
    u = BasisFunction(element)
    v = BasisFunction(element)
    i = Index()
    j = Index()
    w = 0.1*u[i].dx(j)*v[i].dx(j)*dx
    print w
