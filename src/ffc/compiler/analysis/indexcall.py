"""This module provides a simple way of calling functions recursively
on all indices associated with an element of the form algebra"""

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2004-10-13 -- 2007-02-06"
__copyright__ = "Copyright (C) 2004-2007 Anders Logg"
__license__  = "GNU GPL Version 2"

def index_call(object, foo, args = None):
    "Call function foo recursively on all indices"

    # FIXME: Looks like we are missing some objects in the recursion?

    if isinstance(object, Form):
        [index_call(object, foo, args) for m in object.monomials]
    elif isinstance(object, Monomial):
        [index_call(c, foo, args) for c in object.constants]
        [index_call(w, foo, args) for w in object.coefficients]
        [index_call(t, foo, args) for t in object.transforms]
        [index_call(v, foo, args) for v in object.basisfunctions]
    elif isinstance(object, Constant):
        object.number.index_call(foo, args)
    elif isinstance(object, BasisFunction):
        [index_call(i, foo, args) for i in object.component]
        [index_call(d, foo, args) for d in object.derivatives]
    elif isinstance(object, Index):
        foo(object, args)

def __index_add(index, args):
    "Add index to list if index is of given type"
    indices = args[0]
    type = args[1]
    if index.type == type:
        indices += [index]
    return

def __index_add_value(index, args):
    "Add index number to list if index is of given type"
    indices = args[0]
    type = args[1]
    if index.type == type:
        indices += [index.index]
    return

def __index_modify_secondary(index, args):
    "Modify secondary index to index of given type if it appears in the list"
    indices = args[0]
    type = args[1]
    if index.index in indices and index.type == Index.SECONDARY:
        index.type = type
    return

def __index_reassign(index, args):
    "Reassign index from old value to new value"
    iold = args[0]
    inew = args[1]
    type = args[2]
    increment = args[3]
    if index.index == iold and index.type == type:
        index.index = inew
        increment[0] = 1
    return
