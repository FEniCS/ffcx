"""This module provides a simple way of calling functions recursively
on all indices associated with an element of the form algebra"""

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2004-10-13 -- 2007-02-06"
__copyright__ = "Copyright (C) 2004-2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# FFC language modules
import algebra
import tokens
import index

def index_call(object, foo, args = None):
    "Call function foo recursively on all indices"

    if isinstance(object, algebra.Form):
        [index_call(m, foo, args) for m in object.monomials]
    elif isinstance(object, algebra.Monomial):
        [index_call(c, foo, args) for c in object.constants]
        [index_call(w, foo, args) for w in object.coefficients]
        [index_call(t, foo, args) for t in object.transforms]
        [index_call(v, foo, args) for v in object.basisfunctions]
    elif isinstance(object, algebra.Constant):
        index_call(object.number, foo, args)
    elif isinstance(object, tokens.Coefficient):
        index_call(object.n0, foo, args)
        index_call(object.n1, foo, args)
    elif isinstance(object, tokens.Transform):
        index_call(object.index0, foo, args)
        index_call(object.index1, foo, args)
    elif isinstance(object, algebra.BasisFunction):
        # FIXME: Why not call index_call not called on basis function index?
        [index_call(i, foo, args) for i in object.component]
        [index_call(d, foo, args) for d in object.derivatives]
    elif isinstance(object, tokens.Derivative):
        index_call(object.index, foo, args)
    elif isinstance(object, index.Index):
        foo(object, args)
    else:
        raise "RuntimeError", "Don't know how to call index on " + str(object)

def index_add(index, args):
    "Add index to list if index is of given type"
    indices = args[0]
    type = args[1]
    if index.type == type:
        indices += [index]
    return

def index_add_value(index, args):
    "Add index number to list if index is of given type"
    indices = args[0]
    type = args[1]
    if index.type == type:
        indices += [index.index]
    return

def index_modify(index, args):
    "Modify index from one type to another"
    indices = args[0]
    from_type = args[1]
    to_type = args[2]
    if index.index in indices and index.type == from_type:
        index.type = to_type
    return

def index_reassign(index, args):
    "Reassign index from old value to new value"
    from_index = args[0]
    to_index = args[1]
    type = args[2]
    increment = args[3]
    if index.index == from_index and index.type == type:
        index.index = to_index
        increment[0] = 1
    return
