# -*- coding: utf-8 -*-

# Copyright (C) 2017 Garth N. Wells
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

import collections
import hashlib
import os
import pickle

from ufl import FiniteElement, VectorElement, MixedElement
from ufl import interval, triangle, tetrahedron
from ffc.compiler import compile_element


def generate_pybind_code(elements, postfix, includes, pickle_hash):
    """Create a string of C++ code for a factory that returns a
    ufc::finite_element instance.

    """

    switch_header = "  switch (i)\n  {\n"
    switch_entry = "  case {}:\n    return std::make_shared<{}{}{}>();\n"
    switch_end = "  default:\n    \
std::cerr << \"Element not found\" << std::endl; \n    \
return std::shared_ptr<ufc::finite_element>();\n  }"

    # Build switch code block
    switches = ""
    for i, element_name, element_index in elements:
        switches += switch_entry.format(i, element_name, postfix, element_index)

    case_block = switch_header + switches + switch_end

    factory_str = """\
std::shared_ptr<ufc::finite_element> factory(int i)
{{
  {}
}}
""".format(case_block)

    wrapper_header = """\
#include <memory>
#include <pybind11/pybind11.h>
"""

    for include in includes:
        wrapper_header += '#include "' + include + '"\n'

    wrapper_str = """\
{0}
namespace py = pybind11;

{1}

PYBIND11_PLUGIN(factory)
{{
  py::module::import("ufc_wrappers.finite_element");

  // Binding code
  py::module m("factory", "finite element factory");

  // Expose dofmap factory function
  m.def("create_element", &factory, "Create a ufc::finite_element");

  // Add pickle file hasg
  m.def("pickle_hash", [](){{ return \"{2}\"; }},
        "Hash string of corresponding pickle file with Python data");

  return m.ptr();
}}"""

    test = wrapper_str.format(wrapper_header, factory_str, pickle_hash)

    return test


def extract_elements(e, subelement_list=None):
    """Extract sub-elements form a UFL element (includes element itself)

    """
    if not subelement_list:
        subelement_list = []

    # Add element
    subelement_list.append(e)

    # Get subelements
    sub_elements = e.sub_elements()

    for e in sub_elements:
        subelement_list = extract_elements(e, subelement_list)

    return subelement_list


def build_element_tuple(ufl_elements):
    counter = 0
    element_data = collections.defaultdict(list)
    for family, elements in ufl_elements.items():
        for e in elements:
            prefix = "e" + str(hashlib.md5(str(e).encode('utf-8')).hexdigest())
            all_elements = extract_elements(e)
            num_subelements = len(set(all_elements)) - 1
            element_data[family].append((counter, prefix, e, num_subelements))
            counter += 1

    return element_data


def build_ufl_element_list():
    """Build collection of UFL elements for inclusion in factory"""

    lagrange_elements = []
    vector_elements = []
    mixed_elements0 = []
    mixed_elements1 = []
    misc_elements = []

    # Lagrange elements
    for cell in (interval, triangle, tetrahedron):
        for p in range(1, 2):
            lagrange_elements.append(FiniteElement("Lagrange", cell, p))
            lagrange_elements.append(VectorElement("Lagrange", cell, p))
            lagrange_elements.append(FiniteElement("Discontinuous Lagrange", cell, p-1))

    # Vector elements
    for cell in (triangle, tetrahedron):
        for p in range(1, 2):
            vector_elements.append(FiniteElement("RT", cell, p))
            vector_elements.append(FiniteElement("BDM", cell, p))
            vector_elements.append(FiniteElement("N1curl", cell, p))
            vector_elements.append(FiniteElement("N2curl", cell, p))

            vector_elements.append(FiniteElement("Discontinuous Raviart-Thomas", cell, p))

    # Mixed elements
    for cell in (interval, triangle, tetrahedron):
        for p in range(1, 4):
            e0 = FiniteElement("Lagrange", cell, p+1)
            e1 = FiniteElement("Lagrange", cell, p)
            e2 = VectorElement("Lagrange", cell, p+1)

            mixed_elements0.append(MixedElement([e0, e0]))
            mixed_elements0.append(MixedElement([e0, e1]))
            mixed_elements0.append(MixedElement([e1, e0]))
            mixed_elements0.append(MixedElement([e2, e1]))
            mixed_elements0.append(MixedElement([MixedElement([e1, e1]), e0]))

    for cell in (triangle, tetrahedron):
        for p in range(1, 2):
            e0 = FiniteElement("Lagrange", cell, p+1)
            e1 = FiniteElement("Lagrange", cell, p)
            e2 = VectorElement("Lagrange", cell, p+1)
            e3 = FiniteElement("BDM", cell, p)
            e4 = FiniteElement("RT", cell, p+1)
            e5 = FiniteElement("N1curl", cell, p)

            mixed_elements1.append(MixedElement([e1, e2]))
            mixed_elements1.append(MixedElement([e3, MixedElement([e4, e5])]))

    # Misc elements
    for cell in (triangle,):
        for p in range(1, 2):
            misc_elements.append(FiniteElement("HHJ", cell, p))

    for cell in (triangle, tetrahedron):
        misc_elements.append(FiniteElement("CR", cell, 1))
        for p in range(1, 2):
            misc_elements.append(FiniteElement("Regge", cell, p))

    return {"lagrange": lagrange_elements,
            "mixed0": mixed_elements0,
            "mixed1": mixed_elements1,
            "vector":  vector_elements, "misc":  misc_elements}


def compile_elements_to_cpp(elements, outpath):
    """Compile UFL elements to UFC code and write code to file"""

    for family, element_data in elements.items():
        # Build file names and open file
        filename_h = "elements_" + family + ".h"
        filename_cpp = "elements_" + family + ".cpp"
        f_h = open(os.path.join(outpath, filename_h), 'w')
        f_cpp = open(os.path.join(outpath, filename_cpp), 'w')

        for i, prefix, e, num_sub_elements in element_data:
            header, implementation = compile_element(e, prefix=prefix, parameters={"split": True})

            # Adjust include in implementation
            old_include = '#include "' + str(prefix) + '.h\"'
            implementation = implementation.replace(old_include, '#include "' + filename_h + '"', 1)
            f_h.write(header)
            f_cpp.write(implementation)

        # Close file
        f_h.close()
        f_cpp.close()


def build_ufc_wrapper_code(filename_pybind, outpath, elements,
                           pickle_hash):
    """Generate and write to file the pybind11 wrapper code"""

    # Generate C++ code
    compile_elements_to_cpp(elements, outpath)

    # Build C++ file
    element_list = [(i, prefix, e_index) for family, ufl_elements in elements.items() for i, prefix, e, e_index in ufl_elements]

    filenames = ["elements_" + family + ".h" for family, ufl_elements in elements.items()]
    code = generate_pybind_code(element_list, "_finite_element_",
                                filenames, pickle_hash)

    # Current directory of file
    f = open(filename_pybind, 'w')
    f.write(code)
    f.close()


def main():

    # Current directory of file
    dir = os.path.dirname(os.path.realpath(__file__))

    # pybind build directory
    build_dir = os.path.join(dir, "build_pybind")
    if not os.path.exists(build_dir):
        os.makedirs(build_dir)

    # Build list of UFL elements
    ufl_elements = build_ufl_element_list()

    # Build list of element tuples (index, prefix, UFL element,
    # num_unique_subelements)
    elements = build_element_tuple(ufl_elements)

    # Create hash
    elements_hash = hashlib.md5(str(elements).encode('utf-8')).hexdigest()

    # Build and compile wrapper code, and return list of elements
    pathout = os.path.join(dir, build_dir)
    filename_pybind = os.path.join(pathout, 'pybind_wrapper.cpp')
    build_ufc_wrapper_code(filename_pybind, pathout, elements,
                           elements_hash)

    # Save element data to pickle file
    filename = os.path.join(dir, "ffc-test-factory",
                            "element_factory_data.p")
    pickle.dump((elements_hash, elements), open(filename, "wb"))


if __name__ == "__main__":
    main()
