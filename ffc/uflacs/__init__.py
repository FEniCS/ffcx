# -*- coding: utf-8 -*-
# Copyright (C) 2011-2017 Martin Sandve Alnæs
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""This is UFLACS, the UFL Analyser and Compiler System."""

__author__ = u"Martin Sandve Alnæs"

from ffc.uflacs.uflacsrepresentation import compute_integral_ir
from ffc.uflacs.uflacsoptimization import optimize_integral_ir
from ffc.uflacs.uflacsgenerator import generate_integral_code
