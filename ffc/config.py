# Copyright (C) 2019 Garth N. Wells
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import logging
import os
from pathlib import Path

logger = logging.getLogger(__name__)


def get_cache_path(parameters):
    """Get the path for the JIT cache"""

    # Get cache path from parameters, with fallback to current working
    # directory
    cache_dir = parameters.get("cache_dir", Path.cwd() / "compile-cache")

    # Check for cache environment variable
    cache_dir = os.getenv('FENICS_CACHE_DIR', cache_dir)

    return Path(cache_dir).expanduser()
