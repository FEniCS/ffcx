#!/bin/env python
#
# Copyright (C) 2013 Anders Logg
#
# This file is part of FFC.
#
# FFC is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# FFC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with FFC. If not, see <http://www.gnu.org/licenses/>.
#
# First added:  2013-06-24
# Last changed: 2013-06-24
#
# This script finds the latest available reference data tag.

import sys

# Get revisions and references
revs = [rev for rev in sys.argv[1:] if not "ffc-reference-data" in rev]
refs = [ref for ref in sys.argv[1:] if "ffc-reference-data" in ref]

# Pick first available reference
for rev in revs:
    ref = "ffc-reference-data-" + rev
    if ref in refs:
        print ref
        exit(0)

# Unable to find reference data ancestor of HEAD
print "none"
