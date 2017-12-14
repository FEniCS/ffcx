#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""This script is the command-line interface to FFC.

It parses command-line arguments and generates code from input UFL form files.
"""

# Copyright (C) 2004-2017 Anders Logg
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
# Modified by Johan Jansson, 2005.
# Modified by Ola Skavhaug, 2006.
# Modified by Dag Lindbo, 2008.
# Modified by Kristian B. Oelgaard 2010.
# Modified by Martin Sandve Alnæs 2017.

# Python modules.
import sys
import getopt
import cProfile
import re
import string
import os
from os import curdir
from os import path
from os import getcwd

# UFL modules.
from ufl.log import UFLException
from ufl.algorithms import load_ufl_file
import ufl

# FFC modules.
from ffc.log import push_level, pop_level, set_indent, ffc_logger
from ffc.log import DEBUG, INFO, WARNING, ERROR
from ffc.parameters import default_parameters
from ffc import __version__ as FFC_VERSION, get_ufc_signature
from ffc.backends.ufc import __version__ as UFC_VERSION
from ffc.backends.ufc import get_include_path
from ffc.compiler import compile_form, compile_element
from ffc.formatting import write_code
from ffc.errorcontrol import compile_with_error_control


def print_error(msg):
    "Print error message (cannot use log system at top level)."
    print("\n".join(["*** FFC: " + line for line in msg.split("\n")]))


def info_version():
    "Print version number."
    print("""\
This is FFC, the FEniCS Form Compiler, version {0}.
UFC backend version {1}, signature {2}.
For further information, visit https://bitbucket.org/fenics-project/ffc/.

Python {3} on {4}
""".format(FFC_VERSION, UFC_VERSION, get_ufc_signature(),
           sys.version, sys.platform))


def info_usage():
    "Print usage information."
    info_version()
    print("""Usage: ffc [OPTION]... input.form

For information about the FFC command-line interface, refer to
the FFC man page which may invoked by 'man ffc' (if installed).
""")


def compile_ufl_data(ufd, prefix, parameters):
    if parameters["error_control"]:
        code_h, code_c = compile_with_error_control(ufd.forms,
                                                    ufd.object_names,
                                                    ufd.reserved_objects,
                                                    prefix,
                                                    parameters)
    elif len(ufd.forms) > 0:
        code_h, code_c = compile_form(ufd.forms, ufd.object_names,
                                      prefix=prefix,
                                      parameters=parameters)
    else:
        code_h, code_c = compile_element(ufd.elements, prefix=prefix,
                                         parameters=parameters)
    return code_h, code_c


def main(args=None):
    """This is the commandline tool for the python module ffc."""
    if args is None:
        args = sys.argv[1:]

    # Get command-line arguments
    try:
        if "-O" in args:
            args[args.index("-O")] = "-O2"
        opts, args = getopt.getopt(args, "hIVSdvsl:r:f:O:o:q:ep",
                                   ["help", "includes", "version", "signature", "debug", "verbose", "silent",
                                    "language=", "representation=", "optimize=",
                                    "output-directory=", "quadrature-rule=", "error-control",
                                    "profile"])
    except getopt.GetoptError:
        info_usage()
        print_error("Illegal command-line arguments.")
        return 1

    # Check for --help
    if ("-h", "") in opts or ("--help", "") in opts:
        info_usage()
        return 0

    # Check for --includes
    if ("-I", "") in opts or ("--includes", "") in opts:
        print(get_include_path())
        return 0

    # Check for --version
    if ("-V", "") in opts or ("--version", "") in opts:
        info_version()
        return 0

    # Check for --signature
    if ("-S", "") in opts or ("--signature", "") in opts:
        print(get_ufc_signature())
        return 0

    # Check that we get at least one file
    if len(args) == 0:
        print_error("Missing file.")
        return 1

    # Get parameters
    parameters = default_parameters()

    # Choose WARNING as default for script
    parameters["log_level"] = WARNING

    # Set default value (not part of in parameters[])
    enable_profile = False

    # Parse command-line parameters
    for opt, arg in opts:
        if opt in ("-v", "--verbose"):
            parameters["log_level"] = INFO
        elif opt in ("-d", "--debug"):
            parameters["log_level"] = DEBUG
        elif opt in ("-s", "--silent"):
            parameters["log_level"] = ERROR
        elif opt in ("-l", "--language"):
            parameters["format"] = arg
        elif opt in ("-r", "--representation"):
            parameters["representation"] = arg
        elif opt in ("-q", "--quadrature-rule"):
            parameters["quadrature_rule"] = arg
        elif opt == "-f":
            if len(arg.split("=")) == 2:
                (key, value) = arg.split("=")
                default = parameters.get(key)
                if isinstance(default, int):
                    value = int(value)
                elif isinstance(default, float):
                    value = float(value)
                parameters[key] = value
            elif len(arg.split("==")) == 1:
                key = arg.split("=")[0]
                if key.startswith("no-"):
                    key = key[3:]
                    value = False
                else:
                    value = True
                parameters[key] = value
            else:
                info_usage()
                return 1
        elif opt in ("-O", "--optimize"):
            parameters["optimize"] = bool(int(arg))
        elif opt in ("-o", "--output-directory"):
            parameters["output_dir"] = arg
        elif opt in ("-e", "--error-control"):
            parameters["error_control"] = True
        elif opt in ("-p", "--profile"):
            enable_profile = True

    # Set log_level
    push_level(parameters["log_level"])

    # FIXME: This is terrible!
    # Set UFL precision
    #ufl.constantvalue.precision = int(parameters["precision"])

    # Print a versioning message if verbose output was requested
    if parameters["log_level"] <= INFO:
        info_version()

    # Call parser and compiler for each file
    resultcode = 0
    init_indent = ffc_logger._indent_level
    try:
        resultcode = _compile_files(args, parameters, enable_profile)
    finally:
        # Reset logging level and indent
        pop_level()
        set_indent(init_indent)

    return resultcode


def _compile_files(args, parameters, enable_profile):
    # Call parser and compiler for each file
    for filename in args:

        # Get filename prefix and suffix
        prefix, suffix = os.path.splitext(os.path.basename(filename))
        suffix = suffix.replace(os.path.extsep, "")

        # Check file suffix
        if suffix != "ufl":
            print_error("Expecting a UFL form file (.ufl).")
            return 1

        # Remove weird characters (file system allows more than the C
        # preprocessor)
        prefix = re.subn("[^{}]".format(string.ascii_letters + string.digits + "_"), "!", prefix)[0]
        prefix = re.subn("!+", "_", prefix)[0]

        # Turn on profiling
        if enable_profile:
            pr = cProfile.Profile()
            pr.enable()

        # Load UFL file
        ufd = load_ufl_file(filename)

        # Previously wrapped in try-except, disabled to actually get information we need
        #try:

        # Generate code
        code_h, code_c = compile_ufl_data(ufd, prefix, parameters)

        # Write to file
        write_code(code_h, code_c, prefix, parameters)

        #except Exception as exception:
        #    # Catch exceptions only when not in debug mode
        #    if parameters["log_level"] <= DEBUG:
        #        raise
        #    else:
        #        print("")
        #        print_error(str(exception))
        #        print_error("To get more information about this error, rerun FFC with --debug.")
        #        return 1

        # Turn off profiling and write status to file
        if enable_profile:
            pr.disable()
            pfn = "ffc_{0}.profile".format(prefix)
            pr.dump_stats(pfn)
            print("Wrote profiling info to file {0}".format(pfn))

    return 0
