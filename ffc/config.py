# Copyright (C) 2019 Garth N. Wells
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import logging
import os
from pathlib import Path

logger = logging.getLogger(__name__)


def get_cache_path(cache_dir):
    """Get the path for the JIT cache."""
    # Use local cache if cache_dir is None
    if cache_dir is None:
        cache_dir = Path.cwd() / "compile-cache"

    # If cache environment variable is set, override
    cache_dir = os.getenv('FENICS_CACHE_DIR', cache_dir)

    return Path(cache_dir).expanduser()
