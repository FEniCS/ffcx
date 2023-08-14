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
import os
import pathlib
import re
import string

import ufl
from ffcx import __version__ as FFCX_VERSION
from ffcx import compiler
from ffcx.options import FFCX_DEFAULT_OPTIONS, get_options

logger = logging.getLogger("ffcx")

parser = argparse.ArgumentParser(
    description="FEniCS Form Compiler (FFCx, https://fenicsproject.org)"
)
parser.add_argument(
    "--version", action="version", version=f"%(prog)s (version {FFCX_VERSION})"
)
parser.add_argument(
    "-o", "--output-directory", type=str, default=".", help="output directory"
)
parser.add_argument("-L", "--language", type=str, default="C", help="output language")
parser.add_argument("--visualise", action="store_true", help="visualise the IR graph")
parser.add_argument("-p", "--profile", action="store_true", help="enable profiling")

# Add all options from FFCx option system
for opt_name, (opt_val, opt_desc) in FFCX_DEFAULT_OPTIONS.items():
    parser.add_argument(
        f"--{opt_name}", type=type(opt_val), help=f"{opt_desc} (default={opt_val})"
    )

parser.add_argument("ufl_file", nargs="+", help="UFL file(s) to be compiled")


def main(args=None):
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
        prefix = re.subn(
            "[^{}]".format(string.ascii_letters + string.digits + "_"), "!", prefix
        )[0]
        prefix = re.subn("!+", "_", prefix)[0]

        # Turn on profiling
        if xargs.profile:
            pr = cProfile.Profile()
            pr.enable()

        # Load UFL file
        ufd = ufl.algorithms.load_ufl_file(filename)

        # Generate code
        code_header, code_main = compiler.compile_ufl_objects(
            ufd.forms + ufd.expressions + ufd.elements,
            ufd.object_names,
            prefix=prefix,
            options=options,
            visualise=xargs.visualise,
        )

        # Write to file

        suffix_map = {
            "C": (".h", ".c"),
            "c++": (".hpp", ".cpp"),
            "FORTRAN": (None, ".f90"),
            "numba": (None, "_numba.py"),
        }
        head_suffix, main_suffix = suffix_map.get(options["language"], (".h", ".c"))
        output_dir = xargs.output_directory
        if head_suffix:
            filename = os.path.join(output_dir, prefix + head_suffix)
            with open(filename, "w") as fd:
                fd.write(code_header)
        filename = os.path.join(output_dir, prefix + main_suffix)
        with open(filename, "w") as fd:
            fd.write(code_main)

        # Turn off profiling and write status to file
        if xargs.profile:
            pr.disable()
            pfn = f"ffcx_{prefix}.profile"
            pr.dump_stats(pfn)

    return 0
