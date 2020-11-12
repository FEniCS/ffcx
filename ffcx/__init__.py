# Copyright (C) 2009-2018 FEniCS Project
#
# This file is part of FFCX.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

"""FEniCS Form Compiler (FFC).

FFC compiles finite element variational forms into C code.

"""

import logging

import pkg_resources

# Import default parameters
from ffcx.parameters import get_parameters  # noqa: F401

__version__ = pkg_resources.get_distribution("fenics-ffcx").version

logging.basicConfig()
logger = logging.getLogger("ffcx")
logging.captureWarnings(capture=True)
