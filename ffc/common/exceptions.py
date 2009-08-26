"This module defines exceptions for FFC."

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2005-09-14 -- 2009-04-25"
__copyright__ = "Copyright (C) 2005-2009 Anders Logg"
__license__  = "GNU GPL version 3 or any later version"

# FIXME: KBO: This module does not appear to be used anywhere
# find -name "*.py"|xargs grep exceptions
class FFCError(Exception):
    "Base class for FFC exceptions."
    pass
