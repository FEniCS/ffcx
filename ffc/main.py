# -*- coding: utf-8 -*-
# Copyright (C) 2004-2018 Anders Logg and Garth N. Wells
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Command-line interface to FFC.

It parses command-line arguments and generates code from input UFL form files.
"""

import argparse
import cProfile
import logging
import os
import pathlib
import re
import string

import ufl
from ffc import __version__ as FFC_VERSION
from ffc import compiler, formatting
from ffc.parameters import default_parameters

logger = logging.getLogger(__name__)

# Specify command line options
parser = argparse.ArgumentParser(
    description="FEniCS Form Compiler (FFC, https://fenicsproject.org)")
parser.add_argument(
    "-l",
    "--language",
    type=str,
    action='store',
    choices=["ufc", "dolfin"],
    default="ufc",
    help="target language/wrappers")
parser.add_argument(
    "--version",
    action='version',
    version="%(prog)s " + ("(version {})".format(FFC_VERSION)))
parser.add_argument("-d", "--debug", action='store_true', default=False,
                    help="enable debug output")
parser.add_argument("-v", "--verbose", action='store_true', default=False,
                    help="verbose output")
parser.add_argument("-p", "--profile", action='store_true', default=False,
                    help="enable profiling")
parser.add_argument("-o", "--output-directory", type=str, help="output directory")
parser.add_argument(
    "-r",
    "--representation",
    type=str,
    action='store',
    choices=('uflacs', 'tsfc'),
    default="uflacs",
    help="representation ")
parser.add_argument("-q", "--quadrature-rule", default="auto",
                    help="quadrature rule to apply")
parser.add_argument("--quadrature-degree", default="auto",
                    help="quadrature degree to apply")
parser.add_argument('-f', action="append", default=[], dest="f", metavar="parameter=value",
                    help="options passed through to parameter system")
parser.add_argument("ufl_file", nargs='+', help="UFL file(s) to be read from")
# parser.add_argument("--no-evaluate-basis-derivatives", action='store_true',
#                     help="disable generation of code for evaluating basis derivatives")
# parser.add_argument("-O", "--optimize", action='store_true', default=False,
#                    help="apply optimizations during code generation")


def compile_ufl_data(ufd, prefix, parameters):
    if len(ufd.forms) > 0:
        code_h, code_c = compiler.compile_form(
            ufd.forms, ufd.object_names, prefix=prefix, parameters=parameters)
    else:
        code_h, code_c = compiler.compile_element(
            ufd.elements, ufd.object_names, prefix=prefix, parameters=parameters)
    return code_h, code_c


def main(args=None):
    """Commandline tool for FFC."""

    xargs = parser.parse_args(args)

    # Check for --includes
    # if ("-I", "") in opts or ("--includes", "") in opts:
    #     print(ufc.get_include_path())
    #     return 0

    # Check for --signature
    # if ("-S", "") in opts or ("--signature", "") in opts:
    #     print(ufc.get_signature())
    #     return 0

    # Get parameters
    parameters = default_parameters()

    ffc_logger = logging.getLogger("ffc")

    if xargs.debug:
        ffc_logger.setLevel(logging.DEBUG)
    if xargs.verbose:
        ffc_logger.setLevel(logging.INFO)
    parameters["format"] = xargs.language
    parameters["representation"] = xargs.representation
    parameters["quadrature_rule"] = xargs.quadrature_rule
    parameters["quadrature_degree"] = xargs.quadrature_degree
    if xargs.output_directory:
        parameters["output_dir"] = xargs.output_directory
    # parameters["optimize"] = xargs.optimize
    # parameters["no-evaluate_basis_derivatives"] = True
    for p in xargs.f:
        assert len(p.split("=")) == 2
        key, value = p.split("=")
        assert key in parameters
        parameters[key] = value

    # FIXME: This is terrible!
    # Set UFL precision
    # ufl.constantvalue.precision = int(parameters["precision"])

    # Call parser and compiler for each file
    resultcode = _compile_files(xargs.ufl_file, parameters, xargs.profile)
    return resultcode


def _compile_files(args, parameters, enable_profile):
    # Call parser and compiler for each file
    for filename in args:

        file = pathlib.Path(filename)
        # print("TTTT", file.suffix)
        # print("TTTT (1)", file.stem)

        # Get filename and extention string
        prefix, _ = os.path.splitext(os.path.basename(filename))

        # Check file suffix
        if file.suffix != ".ufl":
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
