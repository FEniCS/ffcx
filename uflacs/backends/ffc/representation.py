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

"""The FFC specific backend to the UFLACS form compiler algorithms."""

import numpy

from ufl.algorithms import replace
from ufl.utils.sorting import sorted_by_count

#from uflacs.params import default_parameters
from uflacs.representation.build_uflacs_ir import build_uflacs_ir
from uflacs.elementtables.terminaltables import TableProvider


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
