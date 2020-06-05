# Copyright (C) 2004-2020 Anders Logg, Garth N. Wells and Michal Habera
#
# This file is part of FFCX.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Command-line interface to FFCX.

Parse command-line arguments and generate code from input UFL form files.
"""

import argparse
import cProfile
import logging
import pathlib
import re
import string

import ufl
from ffcx import __version__ as FFCX_VERSION
from ffcx import compiler, formatting
from ffcx.parameters import default_parameters

logger = logging.getLogger("ffcx")

parser = argparse.ArgumentParser(
    description="FEniCS Form Compiler (FFCX, https://fenicsproject.org)")
parser.add_argument(
    "--version", action='version', version="%(prog)s " + ("(version {})".format(FFCX_VERSION)))
parser.add_argument("-v", "--verbosity", action="count", help="verbose output (-vv for more verbosity)")
parser.add_argument("-o", "--output-directory", type=str, default=".", help="output directory")
parser.add_argument("--visualise", action="store_true", help="visualise the IR graph")
parser.add_argument("-p", "--profile", action='store_true', help="enable profiling")

# Add all parameters from FFC parameter system
parameters = default_parameters()
for param_name, param_val in parameters.items():
    parser.add_argument("--{}".format(param_name), default=param_val,
                        type=type(param_val), help="default \"{}\"".format(param_val))

parser.add_argument("ufl_file", nargs='+', help="UFL file(s) to be compiled")


def main(args=None):
    xargs = parser.parse_args(args)

    ffcx_logger = logging.getLogger("ffcx")
    if xargs.verbosity == 1:
        ffcx_logger.setLevel(logging.INFO)
    if xargs.verbosity == 2:
        ffcx_logger.setLevel(logging.DEBUG)

    # Parse all other parameters
    parameters = default_parameters()
    for param_name, param_val in parameters.items():
        parameters[param_name] = xargs.__dict__.get(param_name)

    # Call parser and compiler for each file
    for filename in xargs.ufl_file:
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
        if xargs.profile:
            pr = cProfile.Profile()
            pr.enable()

        # Load UFL file
        ufd = ufl.algorithms.load_ufl_file(filename)

        # Generate code
        if len(ufd.forms) > 0:
            code_h, code_c = compiler.compile_ufl_objects(
                ufd.forms, ufd.object_names, prefix=prefix, parameters=parameters, visualise=xargs.visualise)
        else:
            code_h, code_c = compiler.compile_ufl_objects(
                ufd.elements, ufd.object_names, prefix=prefix, parameters=parameters, visualise=xargs.visualise)

        # Write to file
        formatting.write_code(code_h, code_c, prefix, xargs.output_directory)

        # Turn off profiling and write status to file
        if xargs.profile:
            pr.disable()
            pfn = "ffcx_{0}.profile".format(prefix)
            pr.dump_stats(pfn)

    return 0
