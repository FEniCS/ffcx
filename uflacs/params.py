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

"""Collection of exposed parameters available to tune form compiler algorithms."""


def default_parameters():
    return {
        "enable_profiling": False,
        "enable_factorization": False,  # True, # Fails for hyperelasticity demo in dolfin, needs debugging
        "max_registers": 1024,  # 8 B * 1024 = 8 KB # TODO: Tune this for something very complex
        "score_threshold": 3,  # TODO: Scoring is work in progress and this will change meaning later
    }
