# -*- coding: utf-8 -*-
# Copyright (C) 2011-2017 Martin Sandve Alnæs
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
# along with UFLACS. If not, see <http://www.gnu.org/licenses/>

"""This is UFLACS, the UFL Analyser and Compiler System."""

__author__ = u"Martin Sandve Alnæs"

from ffc.uflacs.uflacsrepresentation import compute_integral_ir
from ffc.uflacs.uflacsoptimization import optimize_integral_ir
from ffc.uflacs.uflacsgenerator import generate_integral_code
