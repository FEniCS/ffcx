__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2004-11-09 -- 2007-01-11"
__copyright__ = "Copyright (C) 2004-2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# FIXME: Move to codegen

class Declaration:

    """A Declaration represents the declaration of a variable,
    including name (type and identifier) and value."""

    def __init__(self, name, value):
        self.name = name
        self.value = value
        self.used = False
        return

    def __str__(self):
        "Pretty print"
        return self.name + " = " + self.value
