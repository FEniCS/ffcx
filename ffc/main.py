# -*- coding: utf-8 -*-
# Copyright (C) 2004-2017 Anders Logg
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""This script is the command-line interface to FFC.

It parses command-line arguments and generates code from input UFL form files.
"""

import cProfile
import getopt
import logging
import os
import re
import string
import sys

import ufl
from ffc import __version__ as FFC_VERSION
from ffc.backends import ufc
from ffc import compiler
from ffc import formatting
from ffc.parameters import default_parameters

logger = logging.getLogger(__name__)


def info_version():
    "Print version number."
    print("""\
This is FFC, the FEniCS Form Compiler, version {0}.
UFC backend version {1}, signature {2}.
For further information, visit https://bitbucket.org/fenics-project/ffc/.

Python {3} on {4}
""".format(FFC_VERSION, ufc.__version__, ufc.get_signature(), sys.version, sys.platform))


def info_usage():
    "Print usage information."
    info_version()
    print("""Usage: ffc [OPTION]... input.form

For information about the FFC command-line interface, refer to
the FFC man page which may invoked by 'man ffc' (if installed).
""")


def compile_ufl_data(ufd, prefix, parameters):
    if len(ufd.forms) > 0:
        code_h, code_c = compiler.compile_form(
            ufd.forms, ufd.object_names, prefix=prefix, parameters=parameters)
    else:
        code_h, code_c = compiler.compile_element(
            ufd.elements, ufd.object_names, prefix=prefix, parameters=parameters)
    return code_h, code_c


def main(args=None):
    """This is the commandline tool for the python module ffc."""
    if args is None:
        args = sys.argv[1:]

    # Get command-line arguments
    print( )
    try:
        if "-O" in args:
            args[args.index("-O")] = "-O2"
        opts, args = getopt.getopt(args, "hIVSdvsl:r:f:O:o:q:ep", [
            "help", "includes", "version", "signature", "debug", "verbose", "silent", "language=",
            "representation=", "optimize=", "output-directory=", "quadrature-rule=",
            "error-control", "profile"
        ])
    except getopt.GetoptError as e:
        info_usage()
        print(e)
        return 1

    # Check for --help
    if ("-h", "") in opts or ("--help", "") in opts:
        info_usage()
        return 0

    # Check for --includes
    if ("-I", "") in opts or ("--includes", "") in opts:
        print(ufc.get_include_path())
        return 0

    # Check for --version
    if ("-V", "") in opts or ("--version", "") in opts:
        info_version()
        return 0

    # Check for --signature
    if ("-S", "") in opts or ("--signature", "") in opts:
        print(ufc.get_signature())
        return 0

    # Check that we get at least one file
    if len(args) == 0:
        logger.error("Missing file.")
        return 1

    # Get parameters
    parameters = default_parameters()

    # Set default value (not part of in parameters[])
    enable_profile = False

    ffc_logger = logging.getLogger("ffc")

    # Parse command-line parameters
    for opt, arg in opts:
        if opt in ("-v", "--verbose"):
            ffc_logger.setLevel(logging.INFO)
        elif opt in ("-d", "--debug"):
            ffc_logger.setLevel(logging.DEBUG)
        elif opt in ("-s", "--silent"):
            ffc_logger.setLevel(logging.ERROR)
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
        elif opt in ("-p", "--profile"):
            enable_profile = True

    # FIXME: This is terrible!
    # Set UFL precision
    # ufl.constantvalue.precision = int(parameters["precision"])

    # Print a versioning message if verbose output was requested
    if logger.getEffectiveLevel() <= logging.INFO:
        info_version()

    # Call parser and compiler for each file
    resultcode = _compile_files(args, parameters, enable_profile)
    return resultcode


def _compile_files(args, parameters, enable_profile):
    # Call parser and compiler for each file
    for filename in args:

        # Get filename prefix and suffix
        prefix, suffix = os.path.splitext(os.path.basename(filename))
        suffix = suffix.replace(os.path.extsep, "")

        # Check file suffix
        if suffix != "ufl":
            logger.error("Expecting a UFL form file (.ufl).")
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
        ufd = ufl.algorithms.load_ufl_file(filename)

        # Previously wrapped in try-except, disabled to actually get information we need
        # try:

        # Generate code
        code_h, code_c = compile_ufl_data(ufd, prefix, parameters)

        # Write to file
        formatting.write_code(code_h, code_c, prefix, parameters)

        # except Exception as exception:
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
