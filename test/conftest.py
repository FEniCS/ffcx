# Copyright (C) 2020 Michal Habera
#
# This file is part of FFCX.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import pytest


@pytest.fixture(scope="module")
def compile_args():
    return ["-O0", "-Wall", "-Werror", "-pedantic-errors"]
