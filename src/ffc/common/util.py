__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2005-02-04 -- 2005-09-07"
__copyright__ = "Copyright (c) 2005 Anders Logg"
__license__  = "GNU GPL Version 2"

def listcopy(l):
    """Create a copy of the list, calling the copy constructor on each
    object in the list (problems when using copy.deepcopy)."""
    if not l:
        return []
    else:
        return [object.__class__(object) for object in l]

def permutations(l):
    "Return a list of all permutations of the given list."
    if len(l) < 1:
        yield l
    for i in range(len(l)):
        pivot = l[i]
        other = l[:i] + l[i+1:]
        for p in permutations(other):
            yield [pivot] + p
    return

if __name__ == "__main__":

    for p in permutations([0, 1, 2, 3]):
        print p
