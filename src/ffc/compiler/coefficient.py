__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2004-11-11"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

# FFC modules
from index import Index

class Coefficient:
    
    """A Coefficient represents the coefficient for a function
    expressed as a linear combination of basis functions.

    A Coefficient holds the following data:

        element - the FiniteElement defining the function space
        number  - a unique index for the represented Function
        index   - summation Index (matching Index of corresponding BasisFunction)"""

    def __init__(self, element, number = None, index = None):
        "Create Coefficient."
        if isinstance(element, Coefficient):
            # Create Coefficient from Coefficient (copy constructor)
            self.element = element.element
            self.number = element.number
            self.index = Index(element.index)
        else:
            # Create Coefficient with given data
            self.element = element
            self.number = number
            self.index = Index(index)
        return

    def __repr__(self):
        "Print nicely formatted representation of Coefficient."
        return "w" + str(self.number) + "_" + str(self.index)

    def indexcall(self, foo, args = None):
        "Call given function on all Indices."
        self.number.indexcall(foo, args)
        self.index.indexcall(foo, args)
        return
