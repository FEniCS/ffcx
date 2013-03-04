# Copyright (C) 2010-2013 Anders Logg, Kristian B. Oelgaard and Marie E. Rognes
#
# This file is part of FFC.
#
# FFC is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# FFC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with FFC. If not, see <http://www.gnu.org/licenses/>.
#
# Modified by Martin Alnaes, 2013
#
# First added:  2010-01-24
# Last changed: 2013-02-14

_test_code = """\
#include "../../ufctest.h"
#include "%s.h"
#include <fstream>

int main(int argc, char * argv[])
{
  const char jsonfilename[] = "%s.json";
  std::ofstream jsonfile(jsonfilename);
  Printer printer(std::cout, jsonfile);
  printer.begin();

%s%s

  printer.end();
  return 0;
}
"""

def generate_test_code(header_file):
    "Generate test code for given header file."

    # Count the number of forms and elements
    prefix = header_file.split(".h")[0]
    generated_code = open(header_file).read()
    num_forms = generated_code.count("class %s_form_" % prefix.lower())
    num_elements = generated_code.count("class %s_finite_element_" % prefix.lower())

    # Generate tests, either based on forms or elements
    if num_forms > 0:
        benchline = "  bool bench = (argc > 1) && argv[1][0] == 'b';\n"
        tests = ['  %s_form_%d f%d; test_form(f%d, bench, %d, printer);' % (prefix.lower(), i, i, i, i)
                 for i in range(num_forms)]
    else:
        benchline = ""
        tests = ['  %s_finite_element_%d e%d; test_finite_element(e%d, %d, printer);' % (prefix.lower(), i, i, i, i)
                 for i in range(num_elements)]

    # Write file
    test_file = open(prefix + ".cpp", "w")
    test_file.write(_test_code % (prefix, prefix, benchline, "\n".join(tests)))
    test_file.close()

