__author__ = "Anders Logg (logg@tti-c.org)"
__version__ = "0.0.1"
__date__ = "2004-09-29"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

from index import Index

class Derivative:

    """A Derivative represents a derivative in either a given fixed
    coordinate direction (in which case it is a tensor of rank 0) or
    one of many possible coordinate directions (in which case it is a
    tensor of rank 1)."""

    def __init__(self, index):
        self.index = Index(index)

    def __repr__(self):
        "Print nicely formatted representation of Derivative."
        return "(d/dx" + str(self.index) + ")"
