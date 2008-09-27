"This module defines exceptions for FFC."

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2005-09-14"
__copyright__ = "Copyright (C) 2005 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

class FFCError(Exception):
    "Base class for FFC exceptions."
    pass

class FormError(FFCError):
    """Exception raised for errors in the definition of a form.
    
    Attributes:
        expression -- input expression in which the error occurred
        message    -- explanation of the error
    """
    def __init__(self, expression, message):
        "Create FormError."
        self.expression = expression
        self.message = message
        return
