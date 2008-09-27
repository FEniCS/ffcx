__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-03-16 -- 2007-03-16"
__copyright__ = "Copyright (C) 2007 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

class Mapping:

    """
    AFFINE is used for H^1 elements,
    CONTRAVARIANT_PIOLA is used for H(div) elements.
    COVARIANT_PIOLA is used for H(curl) elements,
    """
    # Available mapping types
    AFFINE = 0
    CONTRAVARIANT_PIOLA = 1
    COVARIANT_PIOLA = 2 
