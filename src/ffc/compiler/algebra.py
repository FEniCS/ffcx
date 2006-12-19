"""An algebra for multi-linear forms. Objects of the following classes
are elements of the algebra:

    BasisFunction - basic building block
    Product       - product of BasisFunctions
    Sum           - sum of Products
    Function      - linear combination of BasisFunctions
    Constant      - a constant function on the mesh

Each element of the algebra except Constant can be either
scalar or tensor-valued. Elements of the algebra may also
be multi-valued, taking either a ('+') or ('-') value on
interior facets.

The following operations are supported for all elements of the
algebra:

    Binary +      (tensor ranks of operands must match)
    Binary -      (tensor ranks of operands must match)
    Binary *      (both operands must be scalar)
    Unary  +      (operand scalar or tensor-valued)
    Unary  -      (operand scalar or tensor-valued)
    Unary  []     (operand must be tensor-valued)
    Unary  d/dx   (operand scalar or tensor-valued)
    Unary  ()     (operand must be multi-valued, +/-)"""

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2004-09-27 -- 2006-12-19"
__copyright__ = "Copyright (C) 2004-2006 Anders Logg"
__license__  = "GNU GPL Version 2"

# Modified by Garth N. Wells 2006
# Modified by Kristian Oelgaard 2006
# Modified by Marie Rognes 2006

# Python modules
import sys

# FFC common modules
sys.path.append("../../")
from ffc.common.exceptions import *
from ffc.common.debug import *
from ffc.common.util import *

# FFC compiler modules
#from finiteelement import FiniteElement
from reassign import *
from tokens import *
from variables import *
from index import Index
from integral import *
from restriction import *
from piola import *

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
        if isinstance(other, Element):
            return Sum(self) * other
        elif isinstance(other, Integral):
            return Sum(self) * other
        elif isinstance(other, float):
            return Sum(self) * other
        elif isinstance(other, int):
            return Sum(self) * float(other)
        else:
            raise FormError, ((self, other), "Product not defined for given operands.")

    def __rmul__(self, other):
        "Operator: Element * Element (other * self)"
        return Sum(self) * other

    def __div__(self, other):
        "Operator: Element / Element (only works if other is scalar)."
        return Sum(self) * (~other)

    def __rdiv__(self, other):
        "Operator: Element / Element (only works if other is scalar)."
        return Sum(other) * (~self)

    def __invert__(self):
        "Operator: ~Element"
        raise FormError, (self, "Only Constants and numeric constants can be inverted.")

    def __pos__(self):
        "Operator: +Element"
        return Sum(self)

    def __neg__(self):
        "Operator: -Element"
        return -Sum(self)

    def __getitem__(self, component):
        "Operator: Element[component], pick given component."
        return Sum(self)[component]

    def __len__(self):
        "Operator: len(Element)"
        return len(Sum(self))

    def __call__(self, restriction = None):
        "Operator: Element(restriction), restrict multi-valued function."
        return Sum(self)(restriction)

    def dx(self, index = None):
        "Operator: (d/dx)Element in given coordinate direction."
        return Sum(self).dx(index)

    def rank(self):
        "Return value rank of Element."
        return Sum(self).rank()

    def __repr__(self):
        "Print nicely formatted representation of Element."
        return Sum(self).__repr__()

class Constant(Element):
    """A Constant represents a numerical constant or a Function that
    is constant over the mesh.

    Attributes:

        number   - a unique index identifying the Constant.
        inverted - a boolean, true if Constant is inverted
    """

    def __init__(self, constant = None):
        "Create Constant."
        if isinstance(constant, Constant):
            # Create Constant from Constant (copy constructor)
            self.number = Index(constant.number)
            self.inverted = bool(constant.inverted)
        elif constant == None:
            # Create a new number
            self.number = Index("constant")
            self.inverted = False
        else:
            raise FormError, (constant, "Unable to create Constant from given expression.")
        return

    def __invert__(self):
        c = Constant(self)
        if self.inverted:
            c.inverted = False
        else:
            c.inverted = True
        return c

    def __repr__(self):
        "Print nicely formatted representation of Constant."
        return "c" + str(self.number)

    def indexcall(self, foo, args = None):
        "Call given function on all Indices."
        self.number.indexcall(foo, args)
        return
        
class Function(Element):
    """A Function represents a projection of a given function onto a
    finite element space, expressed as a linear combination of
    BasisFunctions.

    Attributes:

        n0 - a unique Index identifying the original function
        n1 - a unique Index identifying the projected function
        e0 - a Finite Element defining the original space
        e1 - a Finite Element defining the projection space
        P  - the projection matrix from e0 to e1

    If projection is not None, then the coefficients of the expansion
    of the function in the current basis should be obtained by applying
    the projection matrix onto the coefficients of the expansion in
    the basis given by the FiniteElement e0.
    """

    def __init__(self, element):
        "Create Function."
        if isinstance(element, Function):
            # Create Function from Function (copy constructor)
            self.n0 = Index(element.n0)
            self.n1 = Index(element.n1)
            self.e0 = element.e0
            self.e1 = element.e1
            self.P  = element.P

        else:
            # Create Function for given FiniteElement
            self.n0 = Index("function")
            self.n1 = Index("projection")
            self.e0 = element
            self.e1 = element
            self.P  = None
        return

    def __repr__(self):
        "Print nicely formatted representation of Function."
        return "w" + str(self.n0)

class BasisFunction(Element):
    """A BasisFunction represents a possibly differentiated component
    of a basis function on the reference cell.

    Attributes:

        element     - a FiniteElement
        index       - a basis Index
        component   - a list of component Indices
        restriction - a flag indicating restriction of a multi-valued function
        derivatives - a list of Derivatives
    """

    def __init__(self, element, index = None):
        "Create BasisFunction."
        if index == None and isinstance(element, BasisFunction):
            # Create BasisFunction from BasisFunction (copy constructor)
            self.element = element.element
            self.index = Index(element.index)
            self.component = listcopy(element.component)
            self.restriction = element.restriction
            self.derivatives = listcopy(element.derivatives)
        elif index == None:
            # Create BasisFunction with primary Index (default)
            self.element = element
            self.index = Index("primary")
            self.component = []
            self.restriction = None
            self.derivatives = []
        else:
            # Create BasisFunction with specified Index
            self.element = element
            self.index = Index(index)
            self.component = []
            self.restriction = None
            self.derivatives = []
        return

    def  __getitem__(self, component):
        "Operator: BasisFunction[component], pick given component."
        if self.element.mapping == "Piola":
            return self.pick_component_piola(component)
        else:
            return self.pick_component_default(component)

    def __len__(self):
        "Operator: len(BasisFunction)"
        if len(self.component) >= self.element.rank():
            raise FormError, (self, "Vector length of scalar expression is undefined.")
        return self.element.tensordim(len(self.component))

    def __call__(self, restriction):
        "Operator: BasisFunction(restriction), restrict multi-valued function."
        if not self.restriction == None:
            raise FormError, ("(" + str(restriction) + ")", "BasisFunction is already restricted.")
        else:
            v = BasisFunction(self)
            if restriction == '+':
                v.restriction = Restriction.PLUS
            elif restriction == '-':
                v.restriction = Restriction.MINUS
        return v

    def __repr__(self):
        "Print nicely formatted representation of BasisFunction."
        d = "".join([d.__repr__() for d in self.derivatives])
        i = self.index.__repr__()

        if self.component:
            c = "[" + ",".join([c.__repr__() for c in self.component]) + "]"
        else:
            c = ""

        if self.restriction == Restriction.PLUS:
            r = "(+)"
        elif self.restriction == Restriction.MINUS:
            r = "(-)"
        else:
            r = ""

        if len(self.derivatives) > 0:
            return "(" + d + "v" + i + c + r + ")"
        else:
            return d + "v" + i + c + r

    def dx(self, index = None):
        "Operator: (d/dx)BasisFunction in given coordinate direction."
        i = Index() # Create new secondary indexF
        w = Product(self)
        w.basisfunctions[0].derivatives.insert(0, Derivative(self.element, i))
        w.transforms.insert(0, Transform(self.element, i, index, self.restriction))
        return w

    def rank(self):
        "Return value rank of BasisFunction."
        if self.component:
            return 0
        else:
            return self.element.rank()

    def eval(self, iindices, aindices, bindices):
        """Evaluate BasisFunction at given indices, returning a tuple consisting
        of (element, number, component, derivative order, derivative indices).
        This tuple uniquely identifies the (possibly differentiated) basis function."""
        vindex = self.index(iindices, aindices, bindices, [])
        cindex = [i(iindices, aindices, bindices, []) for i in self.component]
        dorder = 0
        dindex = [0 for i in range(self.element.shapedim())]
        for d in self.derivatives:
            dindex[d.index(iindices, aindices, bindices, [])] += 1
            dorder += 1
        dindex = tuple(dindex)
        return (self.element, vindex, cindex, dorder, dindex)

    def indexcall(self, foo, args = None):
        "Call given function on all Indices."
        self.index.indexcall(foo, args)
        [i.indexcall(foo, args) for i in self.component]
        [d.indexcall(foo, args) for d in self.derivatives]
        return

    def pick_component_default(self, component):
        "Pick given component of BasisFunction."
        rank = self.element.rank()
        if self.component or rank == 0:
            raise FormError, (self, "Cannot pick component of scalar BasisFunction.")
        w = Product(self)
        if isinstance(component, list):
            if not rank == len(component):
                raise FormError, (component, "Illegal component index, does not match rank.")
            w.basisfunctions[0].component = listcopy(component) 
        else:
            if not rank == 1:
                raise FormError, (component, "Illegal component index, does not match rank.")
            w.basisfunctions[0].component = [Index(component)]        
        return w

    def pick_component_piola(v, component):
        "Pick given component of BasisFunction mapped with the Piola transform."
        rank = v.element.rank()
        if v.component or rank == 0:
            raise FormError, (self, "Cannot pick component of scalar BasisFunction.")    
        w = Product(v)
        j = Index()
        if isinstance(component, list):
            if not rank == len(component):
                raise FormError, (component, "Illegal component index, does not match rank.")
            end = len(component)
            last = component(end)
            w.transforms = [Transform(v.element, j, last, None, -1)] 
            w.basisfunctions[0].component = [component[1:end-1], Index(j)]
            print "The Piola transform is untested in the tensor-valued case."
        else:  
            if not rank == 1:
                raise FormError, (component, "Illegal component index, does not match rank.") 
            w.transforms = [Transform(v.element, j, component, None, -1)] 
            w.basisfunctions[0].component = [Index(j)]    
            w.determinant = -1
        return w

class Product(Element):
    """A Product represents a product of factors, including
    BasisFunctions and Functions.

    Attributes:

        numeric        - a numeric constant (float)
        constants      - a list of Constants
        coefficients   - a list of Coefficients
        transforms     - a list of Transforms
        basisfunctions - a list of BasisFunctions
        integral       - an Integral
        determinant    - the power of the absolute value of 
                         the determinant of the geometry mapping
    """

    def __init__(self, other = None):
        "Create Product."
        if other == None:
            # Create default Product (unity)
            self.numeric = 1.0
            self.constants = []
            self.coefficients = []
            self.transforms = []
            self.basisfunctions = []
            self.determinant = 0
            self.integral = None
        elif isinstance(other, int) or isinstance(other, float):
            # Create Product from scalar
            self.numeric = float(other)
            self.constants = []
            self.coefficients = []
            self.transforms = []
            self.basisfunctions = []
            self.determinant = 0
            self.integral = None
        elif isinstance(other, Function):
            # Create Product from Function
            index = Index()
            self.numeric = 1.0
            self.constants = []
            self.coefficients = [Coefficient(other, index)]
            self.transforms = []
            self.basisfunctions = [BasisFunction(other.e1, index)]
            self.determinant = 0
            self.integral = None
        elif isinstance(other, Constant):
            # Create Product from Constant
            index = Index()
            self.numeric = 1.0
            self.constants = [Constant(other)]
            self.coefficients = []
            self.transforms = []
            self.basisfunctions = []
            self.determinant = 0
            self.integral = None
        elif isinstance(other, BasisFunction):
            # Create Product from BasisFunction
            self.numeric = 1.0
            self.constants = []
            self.coefficients = []
            self.transforms = []
            self.basisfunctions = [BasisFunction(other)]
            self.determinant = 0
            self.integral = None
        elif isinstance(other, Product):
            # Create Product from Product (copy constructor)
            self.numeric = float(other.numeric)
            self.constants = listcopy(other.constants)
            self.coefficients = listcopy(other.coefficients)
            self.transforms = listcopy(other.transforms)
            self.basisfunctions = listcopy(other.basisfunctions)
            self.determinant = other.determinant
            self.integral = other.integral
        else:
            raise FormError, (other, "Unable to create Product from given expression.")

        return
    
    def __mul__(self, other):
        "Operator: Product * Element"
        if isinstance(other, Integral):
            if not self.integral == None:
                raise FormError, (self, "Integrand can only be integrated once.")
            w = Product(self)
            w.integral = Integral(other)
            w.determinant += 1 
            return w
        elif isinstance(other, Sum):
            return Sum(self) * Sum(other)
        else:
            # Create two copies
            w0 = Product(self)
            w1 = Product(other)
            # Check ranks
            if not w0.rank() == w1.rank() == 0:
                raise FormError, (self, "Operands for product must be scalar.")
            # Reassign all complete Indices to avoid collisions
            reassign_complete(w0, Index.SECONDARY)
            reassign_complete(w0, Index.AUXILIARY)
            reassign_complete(w1, Index.SECONDARY)
            reassign_complete(w1, Index.AUXILIARY)
            # Compute product
            w = Product()
            w.numeric = float(w0.numeric * w1.numeric)
            w.constants = listcopy(w0.constants + w1.constants)
            w.coefficients = listcopy(w0.coefficients + w1.coefficients)
            w.transforms = listcopy(w0.transforms + w1.transforms)
            w.basisfunctions = listcopy(w0.basisfunctions + w1.basisfunctions)
            w.determinant = w0.determinant + w1.determinant
            if w0.integral and w1.integral:
                raise FormError, (self, "Integrand can only be integrated once.")
            elif w0.integral:
                w.integral = Integral(w0.integral)
            elif w1.integral:
                w.integral = Integral(w1.integral)
            else:
                w.integral = None;

            return w

    def __neg__(self):
        "Operator: -Product"
        w = Product(self)
        w.numeric = -w.numeric
        return w

    def  __getitem__(self, component):
        "Operator: Product[component], pick given component."
        # Always scalar if product of more than one basis function
        if not len(self.basisfunctions) == 1:
            raise FormError, (self, "Cannot pick component of scalar expression.")
        # Otherwise, return component of first and only BasisFunction
        p = Product(self)
        p.basisfunctions = [] 
        w = Product(self.basisfunctions[0][component]) 
        return w*p

    def __len__(self):
        "Operator: len(Product)"
        # Always scalar if product of more than one basis function
        if not len(self.basisfunctions) == 1:
            raise FormError, (self, "Vector length of scalar expression is undefined.")
        # Otherwise, return length of first and only BasisFunction
        return len(self.basisfunctions[0])

    def __call__(self, r):
        v = Product(self)
        v.basisfunctions = ([w(r) for w in v.basisfunctions])
        # Same restriction for all basis functions so pick first
        restriction = v.basisfunctions[0].restriction
        for i in range(len(v.transforms)):
            v.transforms[i].restriction = restriction
        return v

    def __repr__(self):
        "Print nicely formatted representation of Product."
        if not (self.coefficients or self.transforms or self.basisfunctions):
            return str(self.numeric)
        if self.numeric == -1.0:
            s = "-"
        elif not self.numeric == 1.0:
            s = str(self.numeric)
        else:
            s = ""
        if self.determinant == 0:
            d = ""
        else:
            d = "|det(F)|^" + "(" + str(self.determinant) + ")"
        c = "".join([w.__repr__() for w in self.constants])
        w = "".join([w.__repr__() for w in self.coefficients])
        t = "".join([t.__repr__() for t in self.transforms])
        v = "*".join([v.__repr__() for v in self.basisfunctions])
        if not self.integral == None:
            i = "*" + self.integral.__repr__()
        else:
            i = ""
        return s + c + d + w + t + " | " + v + i

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

    def rank(self):
        "Return value rank of Product."
        if not self.basisfunctions:
            return 0
        if len(self.basisfunctions) > 1:
            for v in self.basisfunctions:
                if not v.rank() == 0:
                    raise FormError, (self, "Illegal rank for BasisFunction of Product (non-scalar).")
        return self.basisfunctions[0].rank()
            
    def indexcall(self, foo, args = None):
        "Call given function on all Indices."
        [c.indexcall(foo, args) for c in self.constants]
        [w.indexcall(foo, args) for w in self.coefficients]
        [t.indexcall(foo, args) for t in self.transforms]
        [v.indexcall(foo, args) for v in self.basisfunctions]
        return

class Sum(Element):
    """A Sum represents a sum of Products. Each Product will be
    compiled separately, since different Products are probably of
    different rank.

    Attributes:

        products - a list of Products (terms) in the Sum
    """
    
    def __init__(self, other = None):
        "Create Sum."
        if other == None:
            # Create default Sum (zero)
            self.products = []
        elif isinstance(other, int) or isinstance(other, float):
            # Create Sum from float
            self.products = [Product(other)]
        elif isinstance(other, BasisFunction):
            # Create Sum from BasisFunction
            self.products = [Product(other)]
        elif isinstance(other, Function):
            # Create Sum from Function
            self.products = [Product(other)]
        elif isinstance(other, Constant):
            # Create Sum from Constant
            self.products = [Product(other)]
        elif isinstance(other, Product):
            # Create Sum from Product
            self.products = [Product(other)]
        elif isinstance(other, Sum):
            # Create Sum from Sum (copy constructor)
            self.products = listcopy(other.products)
        else:
            raise FormError, (other, "Unable to create Sum from given expression.")
        return

    def __add__(self, other):
        "Operator: Sum + Element"
        w0 = Sum(self)
        w1 = Sum(other)
        if not w0.rank() == w1.rank():
            raise FormError, (self, "Operands for addition have non-matching ranks.")
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
        w = Sum()
        w.products = [-p for p in self.products]
        return w

    def  __getitem__(self, component):
        "Operator: indexing, pick given component."
        w = Sum()
        w.products = [p[component] for p in self.products]
        return w

    def __len__(self):
        "Operator: len(Sum)"
        # Check that all terms have the same length
        for i in range(len(self.products) - 1):
            if not len(self.products[i]) == len(self.products[i + 1]):
                raise FormError, (self, "Terms have different vector length.")
        # Return length of first term
        return len(self.products[0])

    def __call__(self, r):
        v = Sum(self)
        v.products = ([w(r) for w in v.products])
        return v

    def __repr__(self):
        "Print nicely formatted representation of Sum."
        return " + ".join([p.__repr__() for p in self.products])

    def dx(self, index = None):
        "Operator: (d/dx)Sum in given coordinate direction."
        w = Sum()
        for p in self.products:
            w = w + p.dx(index)
        return w

    def rank(self):
        "Return value rank of Sum."
        if not self.products:
            return 0
        for j in range(len(self.products) - 1):
            if not self.products[j].rank() == self.products[j + 1].rank():
                raise FormError, (self, "Terms have different rank.")
        return self.products[0].rank()
  
    def indexcall(self, foo, args = None):
        "Call given function on all Indices."
        [p.indexcall(foo, args) for p in self.products]
        return

class TestFunction(BasisFunction):
    """A TestFunction is the BasisFunction with the lowest primary
    index. We simply pick an index lower than all others (-2)."""

    def __init__(self, element):
        index = Index("primary")
        index.index = -2
        BasisFunction.__init__(self, element, index)
        return

class TrialFunction(BasisFunction):
    """A TrialFunction is the BasisFunction with the next lowest primary
    index. We simply pick an index lower than almost all others (-1)."""

    def __init__(self, element):
        index = Index("primary")
        index.index = -1
        BasisFunction.__init__(self, element, index)
        return
