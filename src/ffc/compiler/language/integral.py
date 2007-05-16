__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2004-09-29 -- 2007-03-05"
__copyright__ = "Copyright (C) 2004-2007 Anders Logg"
__license__  = "GNU GPL Version 2"

class Integral:
    """An Integral represents an integral over a mesh entity."""

    # Available integral types
    CELL = 0
    EXTERIOR_FACET = 1
    INTERIOR_FACET = 2

    def __init__(self, type = "cell", sub_domain = 0):
        "Create Integral of given type."
        if isinstance(type, Integral):
            self.type = type.type
            self.sub_domain = type.sub_domain
        elif  type == "cell":
            self.type = self.CELL
            self.sub_domain = sub_domain
        elif type == "exterior facet":
            self.type = self.EXTERIOR_FACET
            self.sub_domain = sub_domain
        elif type == "interior facet":
            self.type = self.INTERIOR_FACET
            self.sub_domain = sub_domain
        else:
            raise RuntimeError, "Unknown integral type " + str(type) + "."
        return

    def __repr__(self):
        "Print nicely formatted representation of Integral."
        if self.type == self.CELL:
            return "dX(%d)" % self.sub_domain
        elif self.type == self.EXTERIOR_FACET:
            return "ds(%d)" % self.sub_domain
        elif self.type == self.INTERIOR_FACET:
            return "dS(%d)" % self.sub_domain

    def __cmp__(self, other):
        "Check if integrals are equal."
        if not isinstance(other, Integral):
            return -1
        if self.type == other.type and self.sub_domain == other.sub_domain:
            return 0
        return -1 # Ignore self > other
        
