__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2005-02-04 -- 2007-04-02"
__copyright__ = "Copyright (C) 2005-2007 Anders Logg"
__license__  = "GNU GPL Version 2"

def pick_first(values):
    "Check that all values are equal and return the value"
    if not values[:-1] == values[1:]:
        raise RuntimeError, "Different values for method %s for objects %s." % (str(method), str(object_list))
    return values[0]

def listcopy(l):
    """Create a copy of the list, calling the copy constructor on each
    object in the list (problems when using copy.deepcopy)."""
    if not l:
        return []
    else:
        return [object.__class__(object) for object in l]

def permutations(l):
    "Return a list of all permutations of the given list"
    if len(l) < 1:
        yield l
    for i in range(len(l)):
        pivot = l[i]
        other = l[:i] + l[i+1:]
        for p in permutations(other):
            yield [pivot] + p
    return

def indent(s, n):
    "Indent each row of the given string s with n spaces"
    indentation = " "*n
    return indentation + ("\n" + indentation).join(s.split("\n"))

if __name__ == "__main__":

    for p in permutations([0, 1, 2, 3]):
        print p

    text = """\
logg@galerkin:~/work/src/fenics/ffc/ffc/src/ffc/common# ls
constants.py   debug.py       exceptions.pyc  progress.py   util.py
constants.pyc  debug.pyc      __init__.py     progress.pyc  util.py.~1.2.~
CVS            exceptions.py  __init__.pyc    #util.py#     util.pyc"""
    print text
    print indent(text, 2)
