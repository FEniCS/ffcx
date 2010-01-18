"""
Code for generating reference data for element
"""

__author__ = "Marie E. Rognes (meg@simula.no)"
__license__  = "GNU GPL version 3 or any later version"

import commands
from ufl import FiniteElement
from ffc.fiatinterface import create_element
from cppcode import element_code, dofmap_code

def run_command(command):
    "Run system command and collect any errors."
    (status, output) = commands.getstatusoutput(command)
    return (status == 0, output)

def generate_reference(element):

    reference = {}

    # Write ufl code to file:
    open("element.ufl", "w").write("element = %s" % element)

    # Compile element with ffc
    c = "ffc element.ufl"
    (ok, output) = run_command(c)
    if not ok:
        print "FFC compilation failed for %s" % element

    # Fetch codes for this element
    fiat_element = create_element(eval(element))
    codes = element_code(fiat_element)

    for (function, code) in codes.iteritems():

        # Create test (c++)
        open("test.cpp", "w").write(code)

        # Compile c++ code
        c = "g++ `pkg-config --cflags ufc-1` -Werror -o test test.cpp"
        (ok, output) = run_command(c)
        if not ok:
            raise Exception, "GCC compilation failed:\n %s" % str(output)

        # Run test
        (ok, result) = run_command("./test")

        # Add result to reference dictionary
        reference[function] = (code, result)

    return reference

elements = ["FiniteElement('CG', 'triangle', 1)"]
references = {}

for element in elements:

    # Generate reference
    references[element] = generate_reference(element)
    print "Reference: ", references[element]


filename = "REFERENCE.py"
file = open(filename, 'w')
file.write("references = %s" % str(references))


