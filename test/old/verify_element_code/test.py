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
    c = "ffc -r tensor element.ufl"
    (ok, output) = run_command(c)
    if not ok:
        print "FFC compilation failed for %s" % element
        failures += [(element, "FFC compilation", output, "")]
        continue

    reference = references[element]

    # For each function
    for function in reference:

        (code, reference_result) = reference[function]

        # TEMPORARY: Changes in generated code...
        code = code.replace("element_0", "element", 1)

        # Write c++ test code to file
        open("test.cpp", "w").write(code)

        # Compile c++ test
        pkg_config = "`pkg-config --cflags ufc-1` `pkg-config --libs dolfin`"
        c = "g++ %s -Werror -o test test.cpp" % pkg_config
        (ok, output) = run_command(c)
        if not ok:
            print "GCC compilation for %s failed:\n %s" \
                   % (function, str(output))
            failures += [(element, "GCC compilation for %s" % function, output, "")]
            continue

        # Run code
        (ok, result) = run_command("./test")
        if not ok:
            print "Running test for for %s failed:\n %s" \
                   % (function, str(result))
            failures += [(element, "Running test for %s" % function, output, "")]
            continue

        # Check result against reference
        if result != reference_result:
            failures += [(element, function, result, reference_result)]
            print
        else:
            print "\tResult matches reference for %s." % (function)


if failures:
    print "\n%d failure(s) [out of possibly %d]:" % (len(failures), len(references.keys())*len(reference.keys()))
    file = open("failures.txt", 'w')
    for (element, function, result, reference_result) in failures:
        msg = "Result does not match reference for %s for %s" % (function, element)
        print "Result does not match reference for %s for %s" % (function, element)
        file.write(msg + ".\n")
        file.write("refere = %s\n" % str(reference_result))
        file.write("result = %s\n\n" % str(result))
    print "\n%d failure(s) [out of possibly %d]:" % (len(failures), len(references.keys())*len(reference.keys()))
    print "Check failures.txt for more info. Good luck and may the force be with you."
    file.close()

else:
    print "All ok!"
