"""This module defines a collection of tokens used to process and
represent multilinear forms, that is, small basic data types used to
build the data structure representing an element of the form algebra."""

__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2004-09-29 -- 2005-11-17"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

# FFC modules
from index import Index

class Coefficient:    
    """A Coefficient represents the coefficient for a function
    expressed as a linear combination of basis functions.

    Attributes:

        n0    - a unique Index identifying the original function
        n1    - a unique Index identifying the projected function
        e0    - a Finite Element defining the original space
        e1    - a Finite Element defining the projection space
        P     - the projection matrix from e0 to e1
        index - an index for summation against corresponding basis function
    """

    def __init__(self, function, index = None):
        "Create Coefficient."
        if isinstance(function, Coefficient):
            # Create Coefficient from Coefficient (copy constructor)
            self.n0 = Index(function.n0)
            self.n1 = Index(function.n1)
            self.e0 = function.e0
            self.e1 = function.e1
            self.P  = function.P
            self.index = Index(function.index)
        else:
            self.n0 = Index(function.n0)
            self.n1 = Index(function.n1)
            self.e0 = function.e0
            self.e1 = function.e1
            self.P  = function.P
            self.index = Index(index)
        return

    def __repr__(self):
        "Print nicely formatted representation of Coefficient."
        return "w" + str(self.n1) + "_" + str(self.index)

    def indexcall(self, foo, args = None):
        "Call given function on all Indices."
        self.n0.indexcall(foo, args)
        self.n1.indexcall(foo, args)
        self.index.indexcall(foo, args)
        return

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

    def indexcall(self, foo, args = None):
        "Call given function on all Indices."
        self.index.indexcall(foo, args)
        return

class Transform:
    """A Transform represents an element of the inverse Jacobian
    matrix of the affine map from the reference cell. With X the
    coordinates on the reference cell mapped to real coordinates x by
    an affine map x = F(X), a Transform represents the partial
    derivative dX/dx."""
    
    def __init__(self, element, index0 = None, index1 = None):
        "Create Transform."
        if isinstance(element, Transform):
            # Create Transform from Transform (copy constructor)
            self.element = element.element
            self.index0 = Index(element.index0)
            self.index1 = Index(element.index1)
        else:
            # Create Transform from given Indices
            self.element = element
            self.index0 = Index(index0)
            self.index1 = Index(index1)
        return

    def __repr__(self):
        "Print nicely formatted representation of Transform."
        return "(dX" + self.index0.__repr__() + "/" + "dx" + self.index1.__repr__() + ")"

    def indexcall(self, foo, args = None):
        "Call given function on all Indices."
        self.index0.indexcall(foo, args)
        self.index1.indexcall(foo, args)
        return

class Integral:
    """An Integral represents an integral over the interior or the
    boundary of the reference cell."""

    def __init__(self, type = "interior"):
        "Create Integral of given type."
        if isinstance(type, Integral):
            # Create Integral from Integral (copy constructor)
            self.type = "" + type.type;
        elif  type == "interior" or type == "boundary":
            # Create Integral of given type
            self.type = "" + type;
        else:
            raise RuntimeError, "Unknown integral type " + str(type) + "."
        return

    def __repr__(self):
        "Print nicely formatted representation of Integral."
        if self.type == "interior":
            return "dX"
        else:
            return "dS"
