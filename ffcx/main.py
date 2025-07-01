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
from collections.abc import Sequence
from typing import Optional

import ufl

from ffcx import __version__ as FFCX_VERSION
from ffcx import compiler, formatting
from ffcx.options import FFCX_DEFAULT_OPTIONS, get_options

logger = logging.getLogger("ffcx")

parser = argparse.ArgumentParser(
    description="FEniCS Form Compiler (FFCx, https://fenicsproject.org)"
)
parser.add_argument("--version", action="version", version=f"%(prog)s (version {FFCX_VERSION})")
parser.add_argument("-o", "--output-directory", type=str, default=".", help="Output directory")
parser.add_argument(
    "-f",
    "--outfile",
    type=str,
    default=None,
    help="Generated code filename stem. Defaults to the stem of the UFL file name.",
)
parser.add_argument(
    "-n",
    "--namespace",
    type=str,
    default=None,
    help="Namespace prefix used in the generated code. Defaults to the stem of the UFL file name.",
)
parser.add_argument("--visualise", action="store_true", help="Visualise the IR graph.")
parser.add_argument("-p", "--profile", action="store_true", help="Enable profiling.")

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


def main(args: Optional[Sequence[str]] = None) -> int:
    """Run ffcx on a UFL file."""
    xargs = parser.parse_args(args)

    # Parse all other options
    priority_options = {k: v for k, v in xargs.__dict__.items() if v is not None}
    options = get_options(priority_options)

    # Call parser and compiler for each file
    for filename in xargs.ufl_file:
        # Remove weird characters (file system allows more than the C
        # preprocessor)
        file_p = pathlib.Path(filename)
        file_p = str(file_p.stem)
        file_p = re.subn("[^{}]".format(string.ascii_letters + string.digits + "_"), "!", file_p)[0]
        file_p = re.subn("!+", "_", file_p)[0]

        if xargs.namespace is None:
            namespace = file_p
        else:
            namespace = xargs.namespace

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
            prefix=namespace,
            visualise=xargs.visualise,
        )

        # Write to file
        if xargs.outfile is None:
            prefix = file_p
        else:
            prefix = xargs.outfile
        formatting.write_code(code_h, code_c, prefix, xargs.output_directory)

        # Turn off profiling and write status to file
        if xargs.profile:
            pr.disable()
            pfn = f"ffcx_{prefix}.profile"
            pr.dump_stats(pfn)

    return 0
