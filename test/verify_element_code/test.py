"""
Tests for generated element and dof map code
"""

__author__ = "Marie E. Rognes (meg@simula.no)"
__license__  = "GNU GPL version 3 or any later version"


# Last modified: 15.01.2010
import commands

from ufl import FiniteElement
from reference_data import reference

from cppcode import *

def run_command(command):
    "Run system command and collect any errors."
    (status, output) = commands.getstatusoutput(command)
    return (status == 0, output)


failures = []
for element in reference:
    print "\nTreating: %s" % element

    # Write ufl code to file:
    open("element.ufl", "w").write("element = %s" % element)

    # Compile element with ffc
    c = "ffc element.ufl"
    (ok, output) = run_command(c)
    if not ok:
        print "FFC compilation failed for %s" % element

    ref = reference[element]

    # For each function
    for function in ref:

        print "Checking %s: " % function,
        code = eval(function)

        # Write test code to file
        open("test.cpp", "w").write(code)

        # Compile file with ufc
        c = "g++ `pkg-config --cflags ufc-1` -Werror -o test test.cpp"
        # c = "g++ `pkg-config --cflags ufc-1` test test.cpp"
        (ok, output) = run_command(c)
        if not ok:
            print output
            raise Exception, "GCC compilation failed"

        # Run code
        (ok, output) = run_command("./test")
        if output != str(ref[function]):
            failures += [(element, function, output)]
            print
        else:
            print output

if failures:
    print "\n%d failure:" % len(failures)
    for (element, function, output) in failures:
        print "%s for %s failed." % (function, element)

else:
    print "All ok!"
