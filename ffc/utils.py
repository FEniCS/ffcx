# -*- coding: utf-8 -*-
# Copyright (C) 2005-2017 Anders Logg
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later


def all_equal(sequence):
    """Check that all items in list are equal."""
    return sequence[:-1] == sequence[1:]
