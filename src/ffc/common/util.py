__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2005-02-04"
__copyright__ = "Copyright (c) 2005 Anders Logg"
__license__  = "GNU GPL Version 2"

def listcopy(list):
    """Create a copy of the list, calling the copy constructor on each
    object in the list (problems when using copy.deepcopy)."""
    if not list:
        return []
    else:
        return [object.__class__(object) for object in list]
