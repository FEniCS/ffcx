# Copyright (C) 2016 Jan Blechta
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

from __future__ import absolute_import

from ffc.log import info
from ffc.representationutils import initialize_integral_ir

from tsfc.driver import compile_integral
import tsfc.kernel_interface.ufc as ufc_interface


def compute_integral_ir(integral_data,
                        form_data,
                        form_id,
                        element_numbers,
                        classnames,
                        parameters):
    "Compute intermediate represention of integral."

    info("Computing tsfc representation")

    # Initialise representation
    ir = initialize_integral_ir("tsfc", integral_data, form_data, form_id)

    # TSFC treats None and unset differently, so remove None values.
    parameters = {k: v for k, v in parameters.items() if v is not None}

    # TSFC has switched to default "coffee" mode
    parameters.setdefault("mode", "vanilla")

    # Store tsfc generated part separately
    ir["tsfc"] = compile_integral(integral_data, form_data, None, parameters,
                                  interface=ufc_interface)

    return ir
