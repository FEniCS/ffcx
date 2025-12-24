# Copyright (C) 2004-2025 Anders Logg, Garth N. Wells and Michal Habera
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Command-line interface to FFCx.

Parse command-line arguments and generate code from input UFL form
files.
"""

import argparse
import cProfile
import logging
import pathlib
import re
import string
from collections.abc import Sequence

import ufl

from ffcx import __version__ as FFCX_VERSION
from ffcx import compiler, formatting
from ffcx.options import FFCX_DEFAULT_OPTIONS, get_options

logger = logging.getLogger("ffcx")

parser = argparse.ArgumentParser(
    description="FEniCS Form Compiler (FFCx, https://fenicsproject.org)"
)
parser.add_argument("--version", action="version", version=f"%(prog)s (version {FFCX_VERSION})")
parser.add_argument("-d", "--dir", type=str, default=".", help="Output directory.")
parser.add_argument(
    "-o",
    "--outfile",
    nargs="*",
    type=str,
    default=None,
    help="Generated code filename stem. Defaults to the stem of the UFL file name.",
)
parser.add_argument(
    "-n",
    "--namespace",
    nargs="*",
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

parser.add_argument(
    "-i",
    "--input",
    nargs="*",
    help="UFL file(s) to be compiled. This option must be used instead "
    "of positional arguments (ufl_file) when using the options -f or -n.",
)
parser.add_argument(
    "ufl_file",
    nargs="*",
    help="UFL file(s) to be compiled. Positional arguments can be used "
    "only when the -f and -n options are not used.",
)


def main(args: Sequence[str] | None = None) -> int:
    """Run ffcx on a UFL file."""
    logging.captureWarnings(capture=True)

    xargs = parser.parse_args(args)

    # Handle UFL files input
    if xargs.input is not None:
        assert len(xargs.ufl_file) == 0, "Unexpected positional arguments with -i option."
        if xargs.namespace is not None:
            assert len(xargs.namespace) == len(xargs.input), (
                "Number of namespaces must match number of input files."
            )
        if xargs.outfile is not None:
            assert len(xargs.outfile) == len(xargs.input), (
                "Number of output files must match number of input files."
            )

        filenames = xargs.input
    else:
        filenames = xargs.ufl_file

    def sanitise_filename(name: str) -> str:
        """Sanitise name by removing non-alphanumeric characters."""
        name_s = pathlib.Path(name).stem
        name_s = re.subn("[^{}]".format(string.ascii_letters + string.digits + "_"), "!", name_s)[0]
        name_s = re.subn("!+", "_", name_s)[0]
        return name_s

    if xargs.namespace is None:
        namespaces = [sanitise_filename(name) for name in filenames]
    else:
        namespaces = xargs.namespace
    if xargs.outfile is None:
        outfiles = [sanitise_filename(name) for name in filenames]
    else:
        outfiles = xargs.outfile

    # Parse all other options
    priority_options = {k: v for k, v in xargs.__dict__.items() if v is not None}
    options = get_options(priority_options)

    for filename, namespace, outfile in zip(filenames, namespaces, outfiles):
        # Turn on profiling
        if xargs.profile:
            pr = cProfile.Profile()
            pr.enable()

        # Load UFL file
        ufd = ufl.algorithms.load_ufl_file(filename)

        # Generate code
        code, suffixes = compiler.compile_ufl_objects(
            ufd.forms + ufd.expressions + ufd.elements,
            options=options,
            object_names=ufd.object_names,
            namespace=namespace,
            visualise=xargs.visualise,
        )

        # Write to file
        formatting.write_code(code, prefix, suffixes, xargs.output_directory)

        # Turn off profiling and write status to file
        if xargs.profile:
            pr.disable()
            pfn = f"ffcx_{namespace}.profile"
            pr.dump_stats(pfn)

    return 0
