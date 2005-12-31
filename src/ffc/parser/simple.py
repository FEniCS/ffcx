"A simple parser for FFC."

__author__ = "Anders Logg (logg@tti-c.org)"
__date__ = "2004-11-15 -- 2005-10-24"
__copyright__ = "Copyright (c) 2004, 2005 Anders Logg"
__license__  = "GNU GPL Version 2"

def parse(filename, language, options):
    "Parse file with given filename, return name of output file."
    print "Parsing " + filename

    # Read input
    infile = open(filename, "r")
    input = infile.read()
    infile.close()

    # Get prefix of file name
    prefix = filename.replace(".form", "")
    
    # Generate output
    output = """\
import sys
sys.path.append("../../")
sys.path.append("./")
from ffc.compiler.compiler import *

name = "%s"

a = None
L = None
M = None
    
%s

compile([a, L, M], name, \"%s\", %s)
""" % (prefix, input, language, options)

    # Write output
    outname = prefix + ".py"
    outfile = open(outname, "w")
    outfile.write(output)
    outfile.close()
    print "Output written to " + outname

    # Return output file name
    return outname
