__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2004-11-11"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

# FFC modules
from index import Index

class Coefficient:
    
    """A Coefficient represents the coefficient for a function
    expressed as a linear combination of basis functions."""

    def __init__(self, element, name = None, index = None):
        "Create Coefficient."
        if isinstance(element, Coefficient):
            # Create Coefficient from Coefficient (copy constructor)
            self.element = element.element
            self.name = element.name
            self.index = Index(element.index)
        else:
            # Create Coefficient with given data
            self.element = element
            self.name = name
            self.index = Index(index)
        return

    def __repr__(self):
        "Print nicely formatted representation of Coefficient."
        return "c" + str(self.index)

    def indexcall(self, foo, args = None):
        "Call given function on all Indices."
        self.index.indexcall(foo, args)
        return
