#!/usr/bin/env python
# Copyright (C) 2017-2017 Martin Sandve Alnæs
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Run ffcx on a UFL file."""

from ffcx.main import main

if __name__ == "__main__":
    import sys

    sys.exit(main())
