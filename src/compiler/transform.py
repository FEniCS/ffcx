__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2004-10-07"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

# FFC modules
from index import Index

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
