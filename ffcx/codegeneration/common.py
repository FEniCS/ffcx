# Copyright (C) 2025 Paul T. Kühner
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

"""Common code for backend implemenations."""

import string


def template_keys(template: str) -> set[str]:
    """Set of expected data keys of a template."""
    return set(fname for _, fname, _, _ in string.Formatter().parse(template) if fname)
