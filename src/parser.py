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
name = "MyPDE"
from form import *
    
dx = Integral("interior")
ds = Integral("boundary")

%s
form = Form(a, name)
form.compile(\"%s\")
""" % (input, language or "")

    # Write output
    outname = filename.rstrip(".form") + ".py"
    outfile = open(outname, "w")
    outfile.write(output)
    outfile.close()
    print "Output generated to " + outname

    # Return output file name
    return outname
