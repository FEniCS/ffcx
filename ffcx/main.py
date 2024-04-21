# Copyright (C) 2004-2020 Anders Logg, Garth N. Wells and Michal Habera
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Command-line interface to FFCx.

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
from ffcx.options import FFCX_DEFAULT_OPTIONS, get_options

logger = logging.getLogger("ffcx")

parser = argparse.ArgumentParser(
    description="FEniCS Form Compiler (FFCx, https://fenicsproject.org)"
)
parser.add_argument("--version", action="version", version=f"%(prog)s (version {FFCX_VERSION})")
parser.add_argument("-o", "--output-directory", type=str, default=".", help="output directory")
parser.add_argument("--visualise", action="store_true", help="visualise the IR graph")
parser.add_argument("-p", "--profile", action="store_true", help="enable profiling")

# Add all options from FFCx option system
for opt_name, (arg_type, opt_val, opt_desc, choices) in FFCX_DEFAULT_OPTIONS.items():
    if isinstance(opt_val, bool):
        parser.add_argument(
            f"--{opt_name}", action="store_true", help=f"{opt_desc} (default={opt_val})"
        )
    else:
        parser.add_argument(
            f"--{opt_name}", type=arg_type, choices=choices, help=f"{opt_desc} (default={opt_val})"
        )

parser.add_argument("ufl_file", nargs="+", help="UFL file(s) to be compiled")


def main(args=None):
    """Run ffcx on a UFL file."""
    xargs = parser.parse_args(args)

    # Parse all other options
    priority_options = {k: v for k, v in xargs.__dict__.items() if v is not None}
    options = get_options(priority_options)

    # Call parser and compiler for each file
    for filename in xargs.ufl_file:
        file = pathlib.Path(filename)

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
        code_h, code_c = compiler.compile_ufl_objects(
            ufd.forms + ufd.expressions + ufd.elements,
            options=options,
            object_names=ufd.object_names,
            prefix=prefix,
            visualise=xargs.visualise,
        )

        # Write to file
        formatting.write_code(code_h, code_c, prefix, xargs.output_directory)

        # Turn off profiling and write status to file
        if xargs.profile:
            pr.disable()
            pfn = f"ffcx_{prefix}.profile"
            pr.dump_stats(pfn)

    return 0
