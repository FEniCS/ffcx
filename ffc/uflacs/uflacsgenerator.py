# Copyright (C) 2013-2014 Martin Alnaes
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

from ffc.log import info
from ffc.representationutils import initialize_integral_code

from ffc.uflacs.backends.ffc.generation import generate_tabulate_tensor_code

def generate_integral_code(ir, prefix, parameters):
    "Generate code for integral from intermediate representation."

    info("Generating code from ffc.uflacs representation")

    # Generate generic ffc code snippets
    code = initialize_integral_code(ir, prefix, parameters)

    # Generate tabulate_tensor body using uflacs algorithms
    uflacs_code = generate_tabulate_tensor_code(ir, prefix, parameters)

    code["tabulate_tensor"] = uflacs_code["tabulate_tensor"]

    # TODO: Use code generation utils here for consistency
    if ir.get("num_cells") is not None:
        code["num_cells"] = "  return %d;" % (ir["num_cells"],)

    code["additional_includes_set"] = set()
    code["additional_includes_set"].update(ir.get("additional_includes_set",()))
    code["additional_includes_set"].update(uflacs_code["additional_includes_set"])

    return code
