__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2004-11-09 -- 2005-10-04"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

class Declaration:

    """A Declaration represents the declaration of a variable,
    including name (type and identifier) and value."""

    def __init__(self, name, value):
        self.name = name
        self.value = value
        self.used = False
        return

    def __repr__(self):
        return self.name + " = " + self.value
