__author__ = "Anders Logg, Kristian B. Oelgaard and Marie E. Rognes"
__date__ = "2010-01-24"
__copyright__ = "Copyright (C) 2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

_test_code = """\
#include "../ufctest.h"
#include "%s.h"

int main()
{
%s

  return 0;
}
"""

def form_test_code(prefix, num_forms):
    "Return code for testing a set of forms."
    tests = ["  %s_form_%d f%d; test_form(f%d);" % (prefix.lower(), i, i, i) for i in range(num_forms)]
    return _test_code % (prefix, "\n".join(tests))

def element_test_code(prefix, num_elements):
    "Return code for testing a set of elements."
    tests = ["%s_element_%d e%d; test_element(e%d);" % (prefix.lower(), i, i, i) for i in range(num_elements)]
    return _test_code % (prefix, "\n".join(tests))
