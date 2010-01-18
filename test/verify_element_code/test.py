"""
Tests for generated element and dof map code
"""

__author__ = "Marie E. Rognes (meg@simula.no)"
__license__  = "GNU GPL version 3 or any later version"


# Last modified: 15.01.2010
import commands

def run_command(command):
    "Run system command and collect any errors."
    (status, output) = commands.getstatusoutput(command)
    return (status == 0, output)

from REFERENCE import references

failures = []
for element in references:
    print "\nTreating: %s" % element

    # Write ufl code to file:
    open("element.ufl", "w").write("element = %s" % element)

    # Compile element with ffc
    c = "ffc element.ufl"
    (ok, output) = run_command(c)
    if not ok:
        print "FFC compilation failed for %s" % element

    reference = references[element]

    # For each function
    for function in reference:

        (code, reference_result) = reference[function]

        # Write test code to file
        open("test.cpp", "w").write(code)

        # Compile file with ufc
        c = "g++ `pkg-config --cflags ufc-1` -Werror -o test test.cpp"
        (ok, output) = run_command(c)
        if not ok:
            raise Exception, "GCC compilation failed"

        # Run code
        (ok, result) = run_command("./test")
        if result != reference_result:
            failures += [(element, function, output)]
            print
        else:
            print "OK: %s: %s" % (function, result)

if failures:
    print "\n%d failure:" % len(failures)
    for (element, function, output) in failures:
        print "%s for %s failed." % (function, element)

else:
    print "All ok!"
