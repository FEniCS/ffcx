"""An algebra for multilinear forms. Objects of the following classes
are elements of the algebra:

    BasisFunction  - basic building block
    Monomial       - monomial product of BasisFunctions
    Form           - sum of Monomials
    Function       - linear combination of BasisFunctions
    Constant       - a constant function on the mesh

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
__date__ = "2004-09-27 -- 2007-03-20"
__copyright__ = "Copyright (C) 2004-2007 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Garth N. Wells 2006
# Modified by Kristian Oelgaard 2006
# Modified by Marie Rognes 2006

# Python modules
import sys
import math

# FFC common modules
from ffc.common.exceptions import *
from ffc.common.debug import *
from ffc.common.utils import *

# FFC FEM modules
from ffc.fem.mapping import *

# FFC language modules
from tokens import *
from index import *
from integral import *
from restriction import *
from reassignment import *

class Element:
    "Base class for elements of the algebra."

    def __add__(self, other):
        "Operator: Element + Element"
        return Form(self) + other

    def __radd__(self, other):
        "Operator: Element + Element (other + self)"
        return Form(self) + other

    def __sub__(self, other):
        "Operator: Element - Element"
        return Form(self) + (-other)

    def __rsub__(self, other):
        "Operator: Element + Element (other - self)"
        return Form(-self) + other

    def __mul__(self, other):
        "Operator: Element * Element"
        if isinstance(other, Element):
            return Form(self) * other
        elif isinstance(other, Integral):
            return Form(self) * other
        elif isinstance(other, float):
            return Form(self) * other
        elif isinstance(other, int):
            return Form(self) * float(other)
        else:
            raise FormError, ((self, other), "Monomial not defined for given operands.")

    def __rmul__(self, other):
        "Operator: Element * Element (other * self)"
        return Form(self) * other

    def __div__(self, other):
        "Operator: Element / Element (only works if other is scalar)."
        return Form(self) * (~other)

    def __rdiv__(self, other):
        "Operator: Element / Element (only works if other is scalar)."
        return Form(other) * (~self)

    def __invert__(self):
        "Operator: ~Element"
        return ~Form(self)

    def __pos__(self):
        "Operator: +Element"
        return Form(self)

    def __neg__(self):
        "Operator: -Element"
        return -Form(self)

    def __getitem__(self, component):
        "Operator: Element[component], pick given component."
        return Form(self)[component]

    def __len__(self):
        "Operator: len(Element)"
        return len(Form(self))

    def __call__(self, restriction = None):
        "Operator: Element(restriction), restrict multi-valued function."
        return Form(self)(restriction)

    def dx(self, index = None):
        "Operator: (d/dx)Element in given coordinate direction."
        return Form(self).dx(index)

    def value_rank(self):
        "Return value rank of Element."
        return Form(self).value_rank()

    def __repr__(self):
        "Print nicely formatted representation of Element."
        return Form(self).__repr__()

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
            self.n0  = Index(element.n0)
            self.n1  = Index(element.n1)
            self.e0  = element.e0
            self.e1  = element.e1
            self.P   = element.P
            self.ops = [op for of in element.ops]

        else:
            # Create Function for given FiniteElement
            self.n0  = Index("function")
            self.n1  = Index("projection")
            self.e0  = element
            self.e1  = element
            self.P   = None
            self.ops = []
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
        if self.element.value_mapping(component) == Mapping.PIOLA:
            return self.pick_component_piola(component)
        else:
            return self.pick_component_default(component)

    def __len__(self):
        "Operator: len(BasisFunction)"
        if len(self.component) >= self.element.value_rank():
            raise FormError, (self, "Vector length of scalar expression is undefined.")
        return self.element.value_dimension(len(self.component))

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
            elif restriction == '+/-':
                v.restriction = Restriction.CONSTANT
            else:
                raise FormError, (self, "Illegal restriction: " + str(restriction))
        return v

    def __repr__(self):
        "Print nicely formatted representation of BasisFunction."
        d = "".join([d.__repr__() for d in self.derivatives])
        i = self.index.__repr__()

        if self.component:
            c = "[" + ", ".join([c.__repr__() for c in self.component]) + "]"
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
        w = Monomial(self)
        w.basisfunctions[0].derivatives.insert(0, Derivative(self.element, i))
        w.transforms.insert(0, Transform(self.element, i, index, self.restriction))
        return w

    def value_rank(self):
        "Return value rank of BasisFunction."
        if self.component:
            return 0
        else:
            return self.element.value_rank()

    def eval(self, iindices, aindices, bindices):
        """Evaluate BasisFunction at given indices, returning a tuple consisting
        of (element, number, component, derivative order, derivative indices).
        This tuple uniquely identifies the (possibly differentiated) basis function."""
        vindex = self.index(iindices, aindices, bindices, [])
        cindex = [i(iindices, aindices, bindices, []) for i in self.component]
        dorder = 0
        dindex = [0 for i in range(self.element.cell_dimension())]
        for d in self.derivatives:
            dindex[d.index(iindices, aindices, bindices, [])] += 1
            dorder += 1
        dindex = tuple(dindex)
        return (self.element, vindex, cindex, dorder, dindex)

    def pick_component_default(self, component):
        "Pick given component of BasisFunction."
        rank = self.element.value_rank()
        if self.component or rank == 0:
            raise FormError, (self, "Cannot pick component of scalar BasisFunction.")
        w = Monomial(self)
        if isinstance(component, list):
            if not rank == len(component):
                raise FormError, (component, "Illegal component index, does not match rank.")
            w.basisfunctions[0].component = listcopy(component) 
        else:
            if not rank == 1:
                raise FormError, (component, "Illegal component index, does not match rank.")
            w.basisfunctions[0].component = [Index(component)]        
        return w

    def pick_component_piola(self, component):
        "Pick given component of BasisFunction mapped with the Piola transform."
        rank = self.element.value_rank()
        if isinstance(component, list):
            if not rank == len(component):
                raise FormError, (component, "Illegal component index, does not match rank.")
            # The Piola transform for the tensor case requires some thought.
            print "The Piola transform is not implemented in the tensor case!"
            return self.pick_component_default(component)
        if not rank == 1:
            raise FormError, (component, "Illegal component index, does not match rank.") 

        (sub_element, offset) = self.element.value_offset(component)
        w = Monomial(self)
        i = Index(component) - offset
        j = Index("secondary", range(self.element.cell_dimension()));
        w.transforms = [Transform(self.element, j, i, self.restriction, Transform.J)] 
        w.basisfunctions[0].component = [j + offset]    
        w.determinants += [Determinant(-1, self.restriction)]
        return w

class Monomial(Element):
    """A Monomial represents a monomial product of factors, including
    BasisFunctions and Functions.

    Attributes:

        numeric        - a numeric constant (float)
        constants      - a list of Constants
        coefficients   - a list of Coefficients
        transforms     - a list of Transforms
        basisfunctions - a list of BasisFunctions
        integral       - an Integral
        determinants    - a list of Determinants
    """

    def __init__(self, other = None):
        "Create Monomial."
        if other == None:
            # Create default Monomial (unity)
            self.numeric = 1.0
            self.constants = []
            self.coefficients = []
            self.transforms = []
            self.basisfunctions = []
            self.determinants = []
            self.integral = None
        elif isinstance(other, int) or isinstance(other, float):
            # Create Monomial from scalar
            self.numeric = float(other)
            self.constants = []
            self.coefficients = []
            self.transforms = []
            self.basisfunctions = []
            self.determinants = []
            self.integral = None
        elif isinstance(other, Function):
            # Create Monomial from Function
            index = Index()
            self.numeric = 1.0
            self.constants = []
            self.coefficients = [Coefficient(other, index)]
            self.transforms = []
            self.basisfunctions = [BasisFunction(other.e1, index)]
            self.determinants = []
            self.integral = None
        elif isinstance(other, BasisFunction):
            # Create Monomial from BasisFunction
            self.numeric = 1.0
            self.constants = []
            self.coefficients = []
            self.transforms = []
            self.basisfunctions = [BasisFunction(other)]
            self.determinants = []
            self.integral = None
        elif isinstance(other, Monomial):
            # Create Monomial from Monomial (copy constructor)
            self.numeric = float(other.numeric)
            self.constants = listcopy(other.constants)
            self.coefficients = listcopy(other.coefficients)
            self.transforms = listcopy(other.transforms)
            self.basisfunctions = listcopy(other.basisfunctions)
            self.determinants = listcopy(other.determinants)
            self.integral = other.integral
        else:
            raise FormError, (other, "Unable to create Monomial from given expression.")

        return
    
    def __mul__(self, other):
        "Operator: Monomial * Element"
        if isinstance(other, Integral):
            if not self.integral == None:
                raise FormError, (self, "Integrand can only be integrated once.")
            w = Monomial(self)
            w.integral = Integral(other)
            return w
        elif isinstance(other, Form):
            return Form(self) * Form(other)
        else:
            # Create two copies
            w0 = Monomial(self)
            w1 = Monomial(other)
            # Check ranks
            if not w0.value_rank() == w1.value_rank() == 0:
                raise FormError, (self, "Operands for monomial must be scalar.")
            # Reassign all complete Indices to avoid collisions
            reassign_complete(w0, Index.SECONDARY)
            reassign_complete(w0, Index.AUXILIARY)
            reassign_complete(w1, Index.SECONDARY)
            reassign_complete(w1, Index.AUXILIARY)
            # Compute monomial
            w = Monomial()
            w.numeric = float(w0.numeric * w1.numeric)
            w.constants = listcopy(w0.constants + w1.constants)
            w.coefficients = listcopy(w0.coefficients + w1.coefficients)
            w.transforms = listcopy(w0.transforms + w1.transforms)
            w.basisfunctions = listcopy(w0.basisfunctions + w1.basisfunctions)
            w.determinants = listcopy(w0.determinants + w1.determinants)
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
        "Operator: -Monomial"
        w = Monomial(self)
        w.numeric = -w.numeric
        return w

    def __invert__(self):
        "Operator: ~Monomial"
        w = Monomial(self)
        w.numeric = 1.0/w.numeric
        for i in range(len(w.coefficients)):
            w.coefficients[i].ops = [Operators.INVERSE] + w.coefficients[i].ops
        return w

    def modulus(self):
        "Take absolute value of monomial"
        w = Monomial(self)
        w.numeric = abs(w.numeric)
        for i in range(len(w.coefficients)):
            w.coefficients[i].ops = [Operators.MODULUS] + w.coefficients[i].ops
        return w

    def sqrt(self):
        "Take square root of monomial"
        w = Monomial(self)
        w.numeric = math.sqrt(w.numeric)
        for i in range(len(w.coefficients)):
            w.coefficients[i].ops = [Operators.SQRT] + w.coefficients[i].ops
        return w

    def  __getitem__(self, component):
        "Operator: Monomial[component], pick given component."
        # Always scalar if monomial of more than one basis function
        if not len(self.basisfunctions) == 1:
            raise FormError, (self, "Cannot pick component of scalar expression.")
        # Otherwise, return component of first and only BasisFunction
        p = Monomial(self)
        p.basisfunctions = [] 
        w = Monomial(self.basisfunctions[0][component]) 
        return w*p

    def __len__(self):
        "Operator: len(Monomial)"
        # Always scalar if monomial of more than one basis function
        if not len(self.basisfunctions) == 1:
            raise FormError, (self, "Vector length of scalar expression is undefined.")
        # Otherwise, return length of first and only BasisFunction
        return len(self.basisfunctions[0])

    def __call__(self, r):
        v = Monomial(self)
        v.basisfunctions = ([w(r) for w in v.basisfunctions])
        # Same restriction for all basis functions so pick first
        restriction = v.basisfunctions[0].restriction
        for i in range(len(v.transforms)):
            v.transforms[i].restriction = restriction
        for i in range(len(v.determinants)):
            v.determinants[i].restriction = restriction
        return v

    def __repr__(self):
        "Print nicely formatted representation of Monomial."
        if not (self.coefficients or self.transforms or self.basisfunctions):
            return str(self.numeric)
        if self.numeric == -1.0:
            s = "-"
        elif not self.numeric == 1.0:
            s = str(self.numeric)
        else:
            s = ""
        d = "".join([d.__repr__() for d in self.determinants])
        c = "".join([w.__repr__() for w in self.constants])
        w = "".join([w.__repr__() for w in self.coefficients])
        t = "".join([t.__repr__() for t in self.transforms])
        v = "*".join([v.__repr__() for v in self.basisfunctions])
        if not self.integral == None:
            i = "*" + str(self.integral)
        else:
            i = ""
        return s + c + d + w + t + " | " + v + i

    def dx(self, index = None):
        "Operator: (d/dx)Monomial in given coordinate direction."
        w = Form()
        for i in range(len(self.basisfunctions)):
            p = Monomial(self)
            p.basisfunctions = []
            for j in range(len(self.basisfunctions)):
                if i == j:
                    p = p * self.basisfunctions[i].dx(index)
                else:
                    p = p * self.basisfunctions[j]
            w = w + p
        return w

    def value_rank(self):
        "Return value rank of Monomial."
        if not self.basisfunctions:
            return 0
        if len(self.basisfunctions) > 1:
            for v in self.basisfunctions:
                if not v.value_rank() == 0:
                    raise FormError, (self, "Illegal rank for BasisFunction of Monomial (non-scalar).")
        return self.basisfunctions[0].value_rank()
            
class Form(Element):
    """A Form represents a sum of Monomials. Each Monomial will be
    compiled separately, since different Monomials are probably of
    different rank.

    Attributes:

        monomials - a list of Monomials (terms) in the Form
    """
    
    def __init__(self, other = None):
        "Create Form."
        if other == None:
            # Create default Form (zero)
            self.monomials = []
        elif isinstance(other, int) or isinstance(other, float):
            # Create Form from float
            self.monomials = [Monomial(other)]
        elif isinstance(other, BasisFunction):
            # Create Form from BasisFunction
            self.monomials = [Monomial(other)]
        elif isinstance(other, Function):
            # Create Form from Function
            self.monomials = [Monomial(other)]
        elif isinstance(other, Monomial):
            # Create Form from Monomial
            self.monomials = [Monomial(other)]
        elif isinstance(other, Form):
            # Create Form from Form (copy constructor)
            self.monomials = listcopy(other.monomials)
        else:
            raise FormError, (other, "Unable to create Form from given expression.")
        return

    def __add__(self, other):
        "Operator: Form + Element"
        w0 = Form(self)
        w1 = Form(other)
        if not w0.value_rank() == w1.value_rank():
            raise FormError, (self, "Operands for addition have non-matching ranks.")
        w = Form()
        w.monomials = w0.monomials + w1.monomials
        return w
    
    def __mul__(self, other):
        "Operator: Form * Element"
        if isinstance(other, Form):
            w = Form()
            w.monomials = [p*q for p in self.monomials for q in other.monomials]
            return w
        else:
            w = Form()
            w.monomials = [p*other for p in self.monomials]
            return w

    def __neg__(self):
        "Operator: -Form"
        w = Form()
        w.monomials = [-p for p in self.monomials]
        return w

    def __invert__(self):
        "Operator: ~Form"
        if len(self.monomials) > 1:
            raise FormError, (self, "Cannot take inverse of sum.")
        w = Form()
        w.monomials = [~p for p in self.monomials]
        return w

    def modulus(self):
        "Take absolute value of form"
        if len(self.monomials) > 1:
            raise FormError, (self, "Cannot take absolute value of sum.")
        w = Form()
        w.monomials = [p.modulus() for p in self.monomials]
        return w

    def sqrt(self):
        "Take square root of form"
        if len(self.monomials) > 1:
            raise FormError, (self, "Cannot take square root of sum.")
        w = Form()
        w.monomials = [p.sqrt() for p in self.monomials]
        return w

    def  __getitem__(self, component):
        "Operator: indexing, pick given component."
        w = Form()
        w.monomials = [p[component] for p in self.monomials]
        return w

    def __len__(self):
        "Operator: len(Form)"
        # Check that all terms have the same length
        for i in range(len(self.monomials) - 1):
            if not len(self.monomials[i]) == len(self.monomials[i + 1]):
                raise FormError, (self, "Terms have different vector length.")
        # Return length of first term
        return len(self.monomials[0])

    def __call__(self, r):
        v = Form(self)
        v.monomials = ([w(r) for w in v.monomials])
        return v

    def __repr__(self):
        "Print nicely formatted representation of Form."
        return " + ".join([p.__repr__() for p in self.monomials])

    def dx(self, index = None):
        "Operator: (d/dx)Form in given coordinate direction."
        w = Form()
        for p in self.monomials:
            w = w + p.dx(index)
        return w

    def value_rank(self):
        "Return value rank of Form."
        if not self.monomials:
            return 0
        for j in range(len(self.monomials) - 1):
            if not self.monomials[j].value_rank() == self.monomials[j + 1].value_rank():
                raise FormError, (self, "Terms have different rank.")
        return self.monomials[0].value_rank()
  
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
