__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2004-11-04"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

class Integral:

    """An Integral represents an integral over the interior or the
    boundary of the reference cell."""

    def __init__(self, type = "interior"):
        "Create Integral of given type."
        if not (type == "interior" or type == "boundary"):
            raise RuntimeError, "Unknown integral type " + str(type) + "."
        self.type = type
        return

    def __repr__(self):
        "Print nicely formatted representation of Integral."
        if self.type == "interior":
            return "dX"
        else:
            return "dS"
