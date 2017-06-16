# -*- coding: utf-8 -*-
"This module defines some basics for generating C++ code."

# Copyright (C) 2009-2017 Anders Logg
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
# Modified by Kristian B. Oelgaard 2011
# Modified by Marie E. Rognes 2010
# Modified by Martin Sandve Alnæs 2013-2017
# Modified by Chris Richardson 2017

# ufc class names

def make_classname(prefix, basename, signature):
    pre = prefix.lower() + "_" if prefix else ""
    sig = str(signature).lower()
    return "%s%s_%s" % (pre, basename, sig)


def make_integral_classname(prefix, integral_type, form_id, subdomain_id):
    basename = "%s_integral_%s" % (integral_type, str(form_id).lower())
    return make_classname(prefix, basename, subdomain_id)
