"""This module defines a collection of tokens used to process and
represent multilinear forms, that is, small basic data types used to
build the data structure representing an element of the form algebra."""

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2004-09-29 -- 2008-08-16"
__copyright__ = "Copyright (C) 2004-2007 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# Modified by Marie Rognes (meg@math.uio.no), 2006

# FFC modules
from index import Index
from restriction import *

class Operators:
    "A collection of simple operators that may be applied to coefficients"

    INVERSE = 0
    MODULUS = 1
    SQRT = 2

class Coefficient:    
    """A Coefficient represents the coefficient for a function
    expressed as a linear combination of basis functions.

    Attributes:

        f     - the function
        n0    - a unique Index identifying the original function
        n1    - a unique Index identifying the projected function
        e0    - a Finite Element defining the original space
        e1    - a Finite Element defining the projection space
        P     - the projection matrix from e0 to e1
        index - an index for summation against corresponding basis function
        ops   - a list of operators to be applied to the function
    """

    def __init__(self, function, index = None):
        "Create Coefficient."
        if isinstance(function, Coefficient):
            # Create Coefficient from Coefficient (copy constructor)
            self.f  = function.f
            self.n0 = Index(function.n0)
            self.n1 = Index(function.n1)
            self.e0 = function.e0
            self.e1 = function.e1
            self.P  = function.P
            self.index = Index(function.index)
            self.ops = [op for op in function.ops]
        else:
            self.f =  function
            self.n0 = Index(function.n0)
            self.n1 = Index(function.n1)
            self.e0 = function.e0
            self.e1 = function.e1
            self.P  = function.P
            self.index = Index(index)
            self.ops = []
        return

    def __repr__(self):
        "Print nicely formatted representation of Coefficient."
        operator_to_string = {Operators.INVERSE: "inv", Operators.MODULUS: "abs", Operators.SQRT: "sqrt"}
        operators = ", ".join([operator_to_string[op] for op in self.ops])
        if not operators == "":
            return "([" + operators + "]" + "w" + str(self.n1) + "_" + str(self.index) + ")"
        else:
            return "w" + str(self.n1) + "_" + str(self.index)

class Derivative:
    """A Derivative represents a derivative on the reference cell in
    either a given fixed coordinate direction (in which case it is a
    tensor of rank 0) or one of many possible coordinate directions
    (in which case it is a tensor of rank 1).

    Attributes:

        element - a FiniteElement
        index   - an Index
    """

    def __init__(self, element, index = None):
        "Create Derivative."
        if isinstance(element, Derivative):
            # Create Derivative from Derivative (copy constructor)
            self.element = element.element
            self.index = Index(element.index)
        else:
            # Create Derivative from given Index
            self.element = element
            self.index = Index(index)
        return

    def __repr__(self):
        "Print nicely formatted representation of Derivative."
        return "(d/dX" + str(self.index) + ")"

class Transform:
    """A Transform represents an element of the (inverse) Jacobian
    matrix of the affine map from the reference cell. With X the
    coordinates on the reference cell mapped to real coordinates x by
    an affine map x = F(X), a Transform represents the partial
    derivative dX/dx (or dx/dX depending on the 'type' of the Transform.)

    Attributes:

        element     - a FiniteElement
        index0      - an Index for variable X
        index1      - an Index for variable x
        type        - the type of the Transform 
        restriction - a Restriction for facet evaluation
    """
    # Available types for the transform
    JINV = 0
    J = 1
    
    def __init__(self, element, index0 = None, index1 = None, restriction = None, type = JINV):
        "Create Transform."
        if isinstance(element, Transform):
            # Create Transform from Transform (copy constructor)
            self.element = element.element
            self.index0 = Index(element.index0)
            self.index1 = Index(element.index1)
            self.type = element.type
            self.restriction = element.restriction
        else:
            # Create Transform from given Indices
            self.element = element
            self.index0 = Index(index0)
            self.index1 = Index(index1)
            self.type = type
            self.restriction = restriction
        return

    def __repr__(self):
        "Print nicely formatted representation of Transform."
        [top, bottom] = ["dX", "dx"]
        [tindex, bindex] = [self.index0.__repr__(), self.index1.__repr__()]
        # If the transform represents dx/dX, we swap notation and indices:
        if self.type == self.J: 
            [bottom, top] = [top, bottom]
            [tindex, bindex] = [bindex, tindex]
            
        restric = ""
        if self.restriction == None:
            restric = ""
        elif self.restriction == Restriction.PLUS:
            restric = "(+)"
        elif self.restriction == Restriction.MINUS:
            restric = "(-)"
        elif self.restriction == Restriction.CONSTANT:
            restric = "(+/-)"
        else:
            raise RuntimeError, "Wrong value for restriction of transform"

        return "(" + top + tindex + "/" + bottom + bindex + ")" + restric

class Determinant:
    """ A Determinant represents a power of the determinant of the
    Jacobean matrix of the affine map from the reference cell.

    Attributes:

        power       - the power of the determinant
        restriction - the restriction of the Jacobean.

    """
    def __init__(self, power, restriction = None):
        "Create Derivative."
        if isinstance(power, Determinant):
            # Create Determinant from Determinant (copy constructor)
            self.power = power.power
            self.restriction = power.restriction
        else:
            # Create Derivative from given power (and restriction)
            self.power = power
            self.restriction = restriction
        return

    def __mul__(self, other):
        if self.restriction == other.restriction:
            return Determinant(self.power + other.power, self.restriction)
        else:
            return None
        
    def __repr__(self):
        "Print nicely formatted representation of Determinant."
        if self.power == 0:
            s = ""
        elif self.power == 1.0:
            s = "det J"
        else:
            s = "(det J)^(%s)" % str(self.power)
        # Maybe make restriction a class and add repr?
        if self.restriction == Restriction.PLUS:
            r = "(+)"
        elif self.restriction == Restriction.MINUS:
            r = "(-)"
        elif self.restriction == Restriction.CONSTANT:
            r = "(+/-)"
        else:
            r = ""
        return s + r
