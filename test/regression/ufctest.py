# -*- coding: utf-8 -*-
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
# Modified by Martin Sandve Alnæs, 2013-2017

_test_code = """\
#include "../../ufctest.h"
#include "{prefix}.h"
#include <fstream>

int main(int argc, char * argv[])
{{
  const char jsonfilename[] = "{prefix}.json";
  std::ofstream jsonfile(jsonfilename);
  Printer printer(jsonfile);
  printer.begin();

{benchline}
{tests}

  printer.end();
  return 0;
}}
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
        tests = ['  {prefix}_form_{i} f{i}; test_form(f{i}, bench, {i}, printer);'.format(prefix=prefix.lower(), i=i)
                 for i in range(num_forms)]
    else:
        benchline = ""
        tests = ['  {prefix}_finite_element_{i} e{i}; test_finite_element(e{i}, {i}, printer);'.format(prefix=prefix.lower(), i=i)
                 for i in range(num_elements)]

    # Write file
    test_file = open(prefix + ".cpp", "w")
    test_file.write(_test_code.format(prefix=prefix, benchline=benchline, tests="\n".join(tests)))
    test_file.close()
