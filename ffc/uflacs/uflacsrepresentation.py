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

from ffc.log import info, error, begin, end, debug_ir, ffc_assert, warning

from ffc.fiatinterface import create_element
from ffc.representationutils import initialize_integral_ir
from ffc.quadrature.parameters import parse_optimise_parameters
from ffc.quadrature.tabulate_basis import tabulate_basis
from ffc.quadrature.quadraturerepresentation import sort_integrals

from ffc.uflacs.backends.ffc.representation import compute_uflacs_integral_ir

def compute_integral_ir(itg_data,
                        form_data,
                        form_id,
                        element_numbers,
                        classnames,
                        parameters):
    "Compute intermediate represention of integral."

    info("Computing uflacs representation")

    # Initialise representation
    ir = initialize_integral_ir("uflacs", itg_data, form_data, form_id)

    # Sort integrals into a dict with quadrature degree and rule as key
    sorted_integrals = sort_integrals(itg_data.integrals,
                                      itg_data.metadata["quadrature_degree"],
                                      itg_data.metadata["quadrature_rule"])

    # TODO: Might want to create the uflacs ir first and then create the tables we need afterwards!
    # Tabulate quadrature points and basis function values in these points
    integrals_dict, psi_tables, quadrature_rules = \
        tabulate_basis(sorted_integrals, form_data, itg_data)

    # Store element numbers, needed for classnames
    ir["element_numbers"] = element_numbers

    # Delegate to flacs to build its intermediate representation and add to ir
    uflacs_ir = compute_uflacs_integral_ir(ir, psi_tables, integrals_dict, form_data, parameters)

    # Store uflacs generated part separately
    ir["uflacs"] = uflacs_ir

    # Create and save the optisation parameters # TODO: Define uflacs specific optimization parameters instead
    #ir["optimise_parameters"] = parse_optimise_parameters(parameters)

    # Save tables for quadrature weights and points
    ir["quadrature_rules"] = quadrature_rules

    # Create dimensions of primary indices, needed to reset the argument 'A'
    # given to tabulate_tensor() by the assembler.
    ir["prim_idims"] = [create_element(ufl_element).space_dimension()
                        for ufl_element in form_data.argument_elements]

    # Added for uflacs, not sure if this is the best way to get this:
    ir["coeff_idims"] = [create_element(ufl_element).space_dimension()
                         for ufl_element in form_data.coefficient_elements]

    return ir
