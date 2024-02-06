# Copyright (C) 2009-2018 FEniCS Project
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""FEniCS Form Compiler (FFCx).

FFCx compiles finite element variational forms into C code.
"""

import importlib.metadata
import logging

# Import default options
from ffcx.options import get_options  # noqa: F401

__version__ = importlib.metadata.version("fenics-ffcx")

logger = logging.getLogger("ffcx")
logging.captureWarnings(capture=True)
