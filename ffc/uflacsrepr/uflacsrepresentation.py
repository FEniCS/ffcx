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
# First added:  2013-02-09
# Last changed: 2013-03-04

#import numpy
from ffc.log import info, error, begin, end, debug_ir, ffc_assert, warning
#from ffc.cpp import format
from ffc.fiatinterface import create_element
from ffc.representationutils import initialize_integral_ir
#from ffc.quadrature.quadratureutils import create_psi_tables
from ffc.quadrature.quadraturerepresentation import _parse_optimise_parameters, _sort_integrals, _tabulate_basis

def compute_integral_ir(itg_data,
                        form_data,
                        form_id,
                        parameters):
    "Compute intermediate represention of integral."

    info("Computing uflacs representation")

    # Initialise representation
    ir = initialize_integral_ir("uflacs", itg_data, form_data, form_id)

    # Sort integrals into a dict with number of integral points as key
    sorted_integrals = _sort_integrals(itg_data.integrals, itg_data.metadata, form_data)

    # Tabulate quadrature points and basis function values in these points
    integrals_dict, psi_tables, quad_weights = \
        _tabulate_basis(sorted_integrals, itg_data.domain_type, form_data)

    # Save tables for quadrature weights and points
    ir["quadrature_weights"] = quad_weights

    # Save tables for basis function values
    ir["psi_tables"] = psi_tables

    # Create dimensions of primary indices, needed to reset the argument 'A'
    # given to tabulate_tensor() by the assembler.
    ir["prim_idims"] = [create_element(ufl_element).space_dimension()
                        for ufl_element in form_data.argument_elements]

    # Added for uflacs, not sure if this is the best way to get this:
    ir["coeff_idims"] = [create_element(ufl_element).space_dimension()
                         for ufl_element in form_data.coefficient_elements]

    # Create and save the optisation parameters.
    ir["optimise_parameters"] = _parse_optimise_parameters(parameters)

    # Delegate to flacs to build its intermediate representation and add to ir
    import uflacs.backends.ffc
    uir = uflacs.backends.ffc.compute_tabulate_tensor_ir(ir, integrals_dict,
                                                         form_data, parameters)
    ir.update(uir)

    return ir
