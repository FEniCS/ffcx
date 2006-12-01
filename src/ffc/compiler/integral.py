__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2004-09-29 -- 2006-12-01"
__copyright__ = "Copyright (C) 2004-2006 Anders Logg"
__license__  = "GNU GPL Version 2"

class Integral:
    """An Integral represents an integral over a mesh entity."""

    # Available integral types
    CELL = 0
    EXTERIOR_FACET = 1
    INTERIOR_FACET = 2

    def __init__(self, type = "cell"):
        "Create Integral of given type."
        if isinstance(type, Integral):
            self.type = type.type
        elif  type == "cell":
            self.type = self.CELL
        elif type == "exterior facet":
            self.type = self.EXTERIOR_FACET
        elif type == "interior facet":
            self.type = self.INTERIOR_FACET
        else:
            raise RuntimeError, "Unknown integral type " + str(type) + "."
        return

    def __repr__(self):
        "Print nicely formatted representation of Integral."
        if self.type == self.CELL:
            return "dX"
        elif self.type == self.EXTERIOR_FACET:
            return "ds"
        elif self.type == self.INTERIOR_FACET:
            return "dS"
