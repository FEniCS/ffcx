#!/usr/bin/env python

# Copyright (C) 2011 Marie E. Rognes
#
# This file is part of UFL.
#
# UFL is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# UFL is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with UFL. If not, see <http://www.gnu.org/licenses/>.
#
# First added:  2011-06-09
# Last changed: 2011-06-09

#
# This is a naive utility script for adding some labels to generated
# .rst and for creating a main level index. It is used by
# scripts/makedoc after generating .rst
#

import os, sys

index_template = """

#############################################
Documentation for FFC v%s
#############################################

FFC Library Reference
=====================

* :ref:`FFC Programmer's Reference <ffc_package>` (Packages and Modules. Everything.)

.. toctree::
   :hidden:
   :maxdepth: 1

   modules


"""

def insert_labels(directory, filenames):
    """
    Insert labels based on filename for those files defined by the
    given filenames relative to directory
    """

    for name in filenames:
        filename = os.path.join(directory, name)
        file = open(filename)
        text = file.read()
        file.close()

        label = "\n.. _%s_package:\n\n" % "_".join(name.split(".")[:-1])
        modded_text = label + text
        print("Adding label to %s" % filename)
        file = open(filename, "w")
        file.write(modded_text)
        file.close()

def generate_index_file(output_dir, version):

    text = index_template % version
    filename = os.path.join(output_dir, "index.rst")

    print("Writing documentation index file to %s" % filename)
    file = open(filename, "w")
    file.write(text)
    file.close()

def main(input_dir, version):

    files = ["ffc.rst", "ffc.tensor.rst", "ffc.quadrature.rst",
             "ffc.errorcontrol.rst"]
    insert_labels(input_dir, files)
    generate_index_file(input_dir, version)

if __name__ == '__main__':

    if len(sys.argv) != 3:
        print("Usage: python generate_index.py input_directory version")
        exit()

    main(sys.argv[1], sys.argv[2])
