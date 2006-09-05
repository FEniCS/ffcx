"ASE output format."

__author__ = "Matthew G. Knepley (knepley@mcs.anl.gov)"
__date__ = "2005-06-03 -- 2006-04-01"
__copyright__ = "Copyright (C) 2005 Matthew G. Knepley"
__license__  = "GNU GPL Version 2"

# Modified by Anders Logg 2005.

try:
    import ASE.Compiler.Python.Cxx
    Cxx = ASE.Compiler.Python.Cxx.Cxx()
    det = Cxx.getVar('jacobianDeterminant')
except ImportError:
    Cxx = None
    det = None

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

# Specify formatting for code generation
format = { "sum": sum,
           "subtract": subtract,
           "multiplication": multiplication,
           "grouping": lambda s: Cxx.getGroup(s),
           "determinant": det,
           "floating point": lambda a: Cxx.getDouble(a),
           "constant": lambda j: None,
           "coefficient": lambda j, k: None,
           "transform": transform,
           "reference tensor" : lambda j, i, a: "(%d, %s, %s)" % (j, str(i), (str(a))),
           "geometry tensor": lambda j, a: 'G%d_%s' % (j, '_'.join([str(index) for index in a])),
           "element tensor": lambda i, k: Cxx.getArrayRef('elementMatrix', k),
           "tmp declaration": lambda j, k: "not defined",
           "tmp access": lambda j, k: "not defined" }

def init(options):
    "Initialize code generation for ase format."
    return

def write(forms, options):
    "Generate code for ase format."
    return
