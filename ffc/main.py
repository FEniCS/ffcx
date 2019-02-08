# -*- coding: utf-8 -*-
# Copyright (C) 2004-2018 Anders Logg and Garth N. Wells
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Command-line interface to FFC.

Parse command-line arguments and generate code from input UFL form files.
"""

import argparse
import cProfile
import logging
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
    "--version", action='version', version="%(prog)s " + ("(version {})".format(FFC_VERSION)))
parser.add_argument("-d", "--debug", action='store_true', help="enable debug output")
parser.add_argument("-v", "--verbose", action='store_true', help="verbose output")
parser.add_argument("-o", "--output-directory", type=str, help="output directory")
parser.add_argument("-p", "--profile", action='store_true', help="enable profiling")
parser.add_argument(
    "-q",
    "--quadrature-rule",
    type=str,
    default="auto",
    help="quadrature rule to apply (default: %(default)s)")
parser.add_argument(
    "--quadrature-degree",
    type=int,
    default=-1,
    help="quadrature degree to apply, auto: -1 (default: %(default)s)")
parser.add_argument(
    "-r",
    "--representation",
    type=str,
    action='store',
    choices=('uflacs', 'tsfc'),
    default="uflacs",
    help="backend to use for compiling forms (default: %(default)s)")
parser.add_argument(
    '-f',
    action="append",
    default=[],
    nargs=2,
    dest="f",
    metavar=("name", "value"),
    help="set existing parameter value in the parameter system, where 'name' is the FFC parameter name")
parser.add_argument(
    '-u',
    action="append",
    default=[],
    nargs=2,
    dest="u",
    metavar=("name", "value"),
    help="add new parameter to the parameter system")
parser.add_argument("ufl_file", nargs='+', help="UFL file(s) to be compiled")


def main(args=None):
    """Commandline tool for FFC."""

    xargs = parser.parse_args(args)
    parameters = default_parameters()
    ffc_logger = logging.getLogger("ffc")

    if xargs.debug:
        ffc_logger.setLevel(logging.DEBUG)
    if xargs.verbose:
        ffc_logger.setLevel(logging.INFO)
    parameters["representation"] = xargs.representation
    parameters["quadrature_rule"] = xargs.quadrature_rule
    parameters["quadrature_degree"] = xargs.quadrature_degree
    if xargs.output_directory:
        parameters["output_dir"] = xargs.output_directory
    for p in xargs.f:
        assert len(p) == 2
        if p[0] not in parameters:
            raise RuntimeError("Command parameter set with -f does not exist in parameters system.")
        parameters[p[0]] = p[1]
    for p in xargs.u:
        assert len(p) == 2
        if p[0] in parameters:
            raise RuntimeError(
                "Command parameter set with -u already exists in parameters system. Use -f.")
        parameters[p[0]] = p[1]

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
        if file.suffix != ".ufl":
            logger.error("Expecting a UFL form file (.ufl).")
            return 1

        # Remove weird characters (file system allows more than the C
        # preprocessor)
        prefix = file.stem
        prefix = re.subn("[^{}]".format(string.ascii_letters + string.digits + "_"), "!", prefix)[0]
        prefix = re.subn("!+", "_", prefix)[0]

        # Turn on profiling
        if enable_profile:
            pr = cProfile.Profile()
            pr.enable()

        # Load UFL file
        ufd = ufl.algorithms.load_ufl_file(filename)

        # Generate code
        if len(ufd.forms) > 0:
            code_h, code_c = compiler.compile_ufl_objects(
                ufd.forms, ufd.object_names, prefix=prefix, parameters=parameters)
        else:
            code_h, code_c = compiler.compile_ufl_objects(
                ufd.elements, ufd.object_names, prefix=prefix, parameters=parameters)

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
