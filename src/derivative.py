__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2004-09-29"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

# FFC modules
from index import Index

class Derivative:
    
    """A Derivative represents a derivative on the reference cell in
    either a given fixed coordinate direction (in which case it is a
    tensor of rank 0) or one of many possible coordinate directions
    (in which case it is a tensor of rank 1)."""

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
