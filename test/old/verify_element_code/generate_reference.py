"""
Code for generating reference data for element
"""

__author__ = "Marie E. Rognes (meg@simula.no)"
__license__  = "GNU GPL version 3 or any later version"

import commands
from ufl import FiniteElement, MixedElement, VectorElement
from cppcode import test_code

def run_command(command):
    "Run system command and collect any errors."
    (status, output) = commands.getstatusoutput(command)
    return (status == 0, output)

def generate_reference(element):
    print "Generating tests for %s" % element
    reference = {}

    # Write ufl code to file:
    open("element.ufl", "w").write("element = %s" % element)

    # Compile element with ffc
    c = "ffc -r tensor element.ufl"
    (ok, output) = run_command(c)
    if not ok:
        print "FFC compilation failed for %s" % element

    # Fetch test code for this element
    try:
        from ffc.fiatinterface import create_element
    except:
        from ffc.createelement import create_element
    codes = test_code(create_element(eval(element)))

    for (function, code) in codes.iteritems():
        # Create test (c++)
        print "\tGenerating, compiling and testing %s" % function
        open("test.cpp", "w").write(code)

        # Compile c++ code
        pkg_config = "`pkg-config --cflags ufc-1` `pkg-config --libs dolfin`"
        c = "g++ %s -Werror -o test test.cpp" % pkg_config
        (ok, output) = run_command(c)
        if not ok:
            raise Exception, "GCC compilation for %s failed:\n %s" \
                  % (function, str(output))

        # Run test
        (ok, result) = run_command("./test")
        if not ok:
            raise Exception, "test for %s did not run!\n" % function

        # Add test code and result to reference dictionary
        reference[function] = (code, result)

    return reference

def look_at_reference():
    from REFERENCE import references

    for (element, reference) in references.iteritems():
        print "-------------------------"
        print element
        for function in reference:
            print "\t %s: %s" % (function, reference[function][1])



if __name__ == "__main__":

    cells = ['triangle', 'tetrahedron']
    families = ["DG", "CG", "RT", "BDM", "N1curl"]
    max_degree = 0
    elements = []

    # List of basic elements to be tested
    elements += ["FiniteElement('%s', '%s', %d)" % (family, cell, k)
                 for family in families
                 for k in range(1, max_degree+1)
                 for cell in cells]
    elements += ["FiniteElement('DG', '%s', 0)" % cell for cell in cells]
    elements += ["FiniteElement('%s', 'interval', %d)" % (family, k)
                 for family in ['DG', 'CG'] for k in range(1, max_degree+1)]

    # Some mixed elements:
    mixed = ["MixedElement([VectorElement('CG', 'triangle', 3), FiniteElement('DG', 'triangle', 2)])",
             "MixedElement([FiniteElement('CG', 'tetrahedron', 2), FiniteElement('DG', 'tetrahedron', 1)])",
             "MixedElement([FiniteElement('DG', 'tetrahedron', 2), FiniteElement('CG', 'tetrahedron', 1)])",
             "MixedElement([FiniteElement('RT', 'tetrahedron', 1), FiniteElement('DG', 'tetrahedron', 0)])",
             "MixedElement([FiniteElement('DG', 'tetrahedron', 1), FiniteElement('BDM', 'tetrahedron', 2)])",
             "MixedElement([FiniteElement('CG', 'triangle', 2), FiniteElement('N1curl', 'triangle', 2)])",
             "MixedElement([FiniteElement('CG', 'tetrahedron', 2), FiniteElement('N1curl', 'tetrahedron', 3)])",
             "MixedElement([VectorElement('CG', 'tetrahedron', 1), FiniteElement('N1curl', 'tetrahedron', 2)])"]

    last = ["MixedElement([VectorElement('CG', 'tetrahedron', 2), FiniteElement('N1curl', 'tetrahedron', 3), FiniteElement('BDM', 'tetrahedron', 3), FiniteElement('DG', 'tetrahedron', 0)])"]

    #elements += mixed

    # Generate reference for each element and store in dictionary
    references = {}
    for element in elements:
        references[element] = generate_reference(element)

    # Dump reference data:
    filename = "REFERENCE.py"
    file = open(filename, 'w')
    file.write("references = %s" % str(references))


