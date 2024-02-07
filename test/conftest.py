# Copyright (C) 2020 Michal Habera
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Test configuration."""

import pytest


@pytest.fixture(scope="module")
def compile_args():
    """Compiler arguments."""
    return ["-O1", "-Wall", "-Werror"]
