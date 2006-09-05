"A simple parser for FFC."

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2004-11-15 -- 2006-03-28"
__copyright__ = "Copyright (C) 2004-2006 Anders Logg"
__license__  = "GNU GPL Version 2"

# FFC common modules
from ffc.common.debug import *

def parse(filename, language, options):
    "Parse file with given filename, return name of output file."
    debug("Parsing " + filename)

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
element = None

%s

if not (a == L == M == None):
  compile([a, L, M], name, \"%s\", %s)
elif not element == None:
  writeFiniteElement(element, name, \"%s\", %s)
else:
  print \"No forms specified, nothing to do.\"

""" % (prefix, input, language, options, language, options)

    # Write output
    outname = prefix + ".py"
    outfile = open(outname, "w")
    outfile.write(output)
    outfile.close()
    debug("Output written to " + outname)

    # Return output file name
    return outname
