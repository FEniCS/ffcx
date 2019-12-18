# Copyright (C) 2009-2018 FEniCS Project
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

"""FEniCS Form Compiler (FFC).

FFC compiles finite element variational forms into C code.

"""

import logging
import pkg_resources

from FIAT import supported_elements

__version__ = pkg_resources.get_distribution("fenics-ffc").version


logging.basicConfig()
logger = logging.getLogger("ffc")
logging.captureWarnings(capture=True)

# Import main function, entry point to script
from ffc.main import main  # noqa: F401

# Import default parameters
from ffc.parameters import default_parameters  # noqa: F401

# Duplicate list of supported elements from FIAT and emove elements from
# list that we don't support or don't trust
supported_elements = sorted(supported_elements.keys())
supported_elements.remove("Argyris")
supported_elements.remove("Hermite")
supported_elements.remove("Morley")
