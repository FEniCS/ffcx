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

import numpy

from ufl.algorithms import replace
from ufl.utils.sorting import sorted_by_count

from ffc.log import info
from ffc.representationutils import initialize_integral_ir
#from ffc.quadrature.parameters import parse_optimise_parameters
from ffc.quadrature.quadraturerepresentation import sort_integrals

# Element table tools
from ffc.fiatinterface import create_element
from ffc.quadrature.tabulate_basis import tabulate_basis
from ffc.uflacs.elementtables.terminaltables import TableProvider

#from ffc.uflacs.params import default_parameters
from ffc.uflacs.representation.build_uflacs_ir import build_uflacs_ir


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

    # Store element numbers TODO: Still in use?
    ir["element_numbers"] = element_numbers

    # Store element classnames
    ir["classnames"] = classnames
    
    # Delegate to flacs to build its intermediate representation and add to ir
    uflacs_ir = compute_uflacs_integral_ir(ir, psi_tables, integrals_dict, form_data, parameters)

    # Store uflacs generated part separately
    ir["uflacs"] = uflacs_ir

    # Create and save the optimisation parameters
    # TODO: Define uflacs specific optimization parameters instead
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

    # TODO: Can easily add more element data here or move above if needed in uflacs
    #unique_elements = element_numbers.keys()
    #ir["fiat_elements"] = { ufl_element: create_element(ufl_element)
    #                        for ufl_element in unique_elements }
    #ir["element_dimensions"] = { ufl_element: fiat_element.space_dimension()
    #                             for ufl_element, fiat_element in ir["fiat_elements"].items() }

    return ir


def compute_uflacs_integral_ir(ir, psi_tables,
                               integrals_dict, form_data,
                               parameters):
    # TODO: Hack before we get default parameters properly into ffc
    #p = default_parameters()
    #p.update(parameters)
    #parameters = p

    # Build coefficient numbering for UFC interface here, to avoid
    # renumbering in UFL and application of replace mapping
    if True:
        # Using the mapped coefficients, numbered by UFL
        coefficient_numbering = {}
        sorted_coefficients = sorted_by_count(form_data.function_replace_map.keys())
        for i, f in enumerate(sorted_coefficients):
            g = form_data.function_replace_map[f]
            coefficient_numbering[g] = i
            assert i == g.count()

        # Replace coefficients so they all have proper element and domain for what's to come
        # TODO: We can avoid the replace call when proper Expression support is in place
        #       and element/domain assignment is removed from compute_form_data.
        integrands = {
            num_points: replace(integrals_dict[num_points].integrand(), form_data.function_replace_map)
            for num_points in sorted(integrals_dict.keys())
            }
    else:
        pass
        #coefficient_numbering = {}
        #coefficient_element = {}
        #coefficient_domain = {}
        #sorted_coefficients = sorted_by_count(form_data.function_replace_map.keys())
        #for i, f in enumerate(sorted_coefficients):
        #    g = form_data.function_replace_map[f]
        #    coefficient_numbering[f] = i
        #    coefficient_element[f] = g.ufl_element()
        #    coefficient_domain[f] = g.ufl_domain()
        #integrands = {
        #    num_points: integrals_dict[num_points].integrand()
        #    for num_points in sorted(integrals_dict.keys())
        #    }
        # then pass coefficient_element and coefficient_domain to the uflacs ir as well

    # Hiding ffc data behind interface that we can
    # improve later to build tables on the fly instead of
    # precomputing psi_tables in ffc, somewhat disconnecting the
    # uflacs representation building from the psi_tables format
    table_provider = TableProvider(psi_tables, parameters)

    # Some more form_data info that we may need to
    # insert in the uflacs_ir but currently don't use
    #form_data.name
    #form_data.coefficient_names
    #form_data.argument_names
    #form_data.integration_domains[0].ufl_cell()
    #form_data.function_replace_map

    # Build the uflacs-specific intermediate representation
    return build_uflacs_ir(ir, integrands, coefficient_numbering, table_provider)
