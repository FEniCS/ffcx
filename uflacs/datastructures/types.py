# -*- coding: utf-8 -*-
# Copyright (C) 2011-2015 Martin Sandve Alnæs
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

"""Basic data structures."""

import numpy

def sufficient_int_type(maxvalue):
    if maxvalue < 2 ** 7:
        dtype = numpy.int8
    elif maxvalue < 2 ** 15:
        dtype = numpy.int16
    elif maxvalue < 2 ** 31:
        dtype = numpy.int32
    else:
        dtype = numpy.int64
    return dtype

def sufficient_uint_type(maxvalue):
    if maxvalue < 2 ** 8:
        dtype = numpy.uint8
    elif maxvalue < 2 ** 16:
        dtype = numpy.uint16
    elif maxvalue < 2 ** 32:
        dtype = numpy.uint32
    else:
        dtype = numpy.uint64
    return dtype
