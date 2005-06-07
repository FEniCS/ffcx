"ASE output format."

__author__ = "Matthew G. Knepley (knepley@mcs.anl.gov)"
__date__ = "2005-06-03"
__copyright__ = "Copyright (c) 2005 Matthew G. Knepley"
__license__  = "GNU GPL Version 2"

try:
    import ASE.Compiler.Python.Cxx
    Cxx = ASE.Compiler.Python.Cxx.Cxx()
except ImportError:
    Cxx = None

dim = 2

def sum(l):
    if not l:
        return None
    import ASE.Compiler.Cxx.Addition
    add = ASE.Compiler.Cxx.Addition.Addition()
    add.setChildren(map(Cxx.getValue, l))
    return add

def subtract(l):
    if not l:
        return None
    import ASE.Compiler.Cxx.Subtraction
    sub = ASE.Compiler.Cxx.Subtraction.Subtraction()
    sub.setChildren(map(Cxx.getValue, l))
    return sub

def multiplication(l):
    if not l:
        return None
    import ASE.Compiler.Cxx.Multiplication
    mult = ASE.Compiler.Cxx.Multiplication.Multiplication()
    mult.setChildren(map(Cxx.getValue, l))
    return mult

def transform(j, k):
    return Cxx.getArrayRef('JInv', int(j)*dim + int(k))

format = { "sum": sum,
           "subtract": subtract,
           "multiplication": multiplication,
           "grouping": lambda s: Cxx.getGroup(s),
           "determinant": Cxx.getVar('jacobianDeterminant'),
           "floating point": lambda a: Cxx.getDouble(a),
           "constant": lambda j: "not defined",
           "coefficient": lambda j, k: "not defined",
           "transform": transform,
           "reference tensor" : lambda j, i, a: "(%d, %s, %s)" % (j, str(i), (str(a))),
           "geometry tensor": lambda j, a: 'G%d_%s' % (j, '_'.join([str(index) for index in a])),
           "element tensor": lambda i, k: Cxx.getArrayRef('elementMatrix', k)}

def compile(forms, license):
    "Generate code for raw format."
    return
