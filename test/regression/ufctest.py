__author__ = "Anders Logg, Kristian B. Oelgaard and Marie E. Rognes"
__date__ = "2010-01-24"
__copyright__ = "Copyright (C) 2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Last changed: 2010-05-11

_test_code = """\
#include "../ufctest.h"
#include "%s.h"

int main()
{
%s

  return 0;
}
"""

def generate_test_code(header_file, bench):
    "Generate test code for given header file."

    # Count the number of forms and elements
    prefix = header_file.split(".h")[0]
    generated_code = open(header_file).read()
    num_forms = generated_code.count("class %s_form_" % prefix.lower())
    num_elements = generated_code.count("class %s_finite_element_" % prefix.lower())

    # Generate tests, either based on forms or elements
    if num_forms > 0:
        tests = ["  %s_form_%d f%d; test_form(f%d, %d);" % (prefix.lower(), i, i, i, bench) for i in range(num_forms)]
    else:
        tests = ["  %s_finite_element_%d e%d; test_finite_element(e%d);" % (prefix.lower(), i, i, i) for i in range(num_elements)]

    # Write file
    test_file = open(prefix + ".cpp", "w")
    test_file.write(_test_code % (prefix, "\n".join(tests)))
    test_file.close()
