"A simple parser for FFC."

__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2004-11-15"
__copyright__ = "Copyright (c) 2004 Anders Logg"
__license__  = "GNU GPL Version 2"

def parse(filename, language):
    "Parse file with given filename, return name of output file."
    print "Parsing " + filename

    # Read input
    infile = open(filename, "r")
    input = infile.read()
    infile.close()
    
    # Generate output
    output = """\
import sys
sys.path.append("../")
sys.path.append("../../")
from ffc.compiler.compiler import *

name = "MyPDE"

a = None
L = None
    
dx = Integral("interior")
ds = Integral("boundary")

i = Index()
j = Index()
k = Index()
l = Index()
m = Index()
n = Index()

%s
compile([a, L], name, \"%s\")
""" % (input, language or "")

    # Write output
    # FIXME: Why does this remove the last 'm' in poissonsystem.form?
    outname = filename.rstrip(".form") + ".py"
    outfile = open(outname, "w")
    outfile.write(output)
    outfile.close()
    print "Output written to " + outname

    # Return output file name
    return outname
