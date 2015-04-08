# Copyright (C) 2011-2015 Martin Sandve Alnes
#
# This file is part of UFLACS.
#
# UFLACS is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# UFLACS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with UFLACS. If not, see <http://www.gnu.org/licenses/>.

"""Precedence list."""

import ufl


def build_precedence_list():
    "Builds a list of operator types by precedence order in the C language."
    # FIXME: Add all types we need here.
    pl = []
    pl.append((ufl.classes.Conditional,))
    pl.append((ufl.classes.OrCondition,))
    pl.append((ufl.classes.AndCondition,))
    pl.append((ufl.classes.EQ, ufl.classes.NE))
    pl.append((ufl.classes.Condition,))  # <,>,<=,>=
    pl.append((ufl.classes.NotCondition,))  # FIXME
    pl.append((ufl.classes.Sum,))
    pl.append((ufl.classes.Product, ufl.classes.Division,))
    # The highest precedence items will never need
    # parentheses around them or their operands
    pl.append((ufl.classes.Power, ufl.classes.MathFunction, ufl.classes.Abs, ufl.classes.BesselFunction,
               ufl.classes.Indexed, ufl.classes.Grad,
               ufl.classes.PositiveRestricted, ufl.classes.NegativeRestricted,
               ufl.classes.Terminal))
    # FIXME: Write a unit test that checks this list against all ufl classes
    return pl


def build_precedence_map():
    from ufl.precedence import build_precedence_mapping
    pm, missing = build_precedence_mapping(build_precedence_list())
    if 0 and missing:  # Enable to see which types we are missing
        print("Missing precedence levels for the types:")
        print("\n".join('  %s' % c for c in missing))
    return pm
