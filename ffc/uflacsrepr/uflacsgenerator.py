# Copyright (C) 2013 Martin Alnaes
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
# First added:  2013-02-12
# Last changed: 2013-02-14

from ffc.representationutils import initialize_integral_code
from ffc.log import info, error, begin, end, debug_ir, ffc_assert, warning
from ffc.cpp import format

def generate_integral_code(ir, prefix, parameters):
    "Generate code for integral from intermediate representation."

    info("Generating code from uflacs representation")

    # Generate generic ffc code snippets
    code = initialize_integral_code(ir, prefix, parameters)

    # Reusing some functions from the quadrature representation code generation
    import ffc.quadrature.quadraturegenerator as qr

    # FIXME: Generate code for basis function tables

    # Delegate to uflacs to generate tabulate_tensor body
    import uflacs.backends.ffc
    ttcode = uflacs.backends.ffc.generate_tabulate_tensor_code(ir, parameters)

    code["tabulate_tensor"] = ttcode

    code["tabulate_tensor_quadrature"] = format["do nothing"] # TODO: Remove
    return code
