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
# Last changed: 2013-02-12

from ffc.representationutils import initialize_integral_ir
from ffc.log import info, error, begin, end, debug_ir, ffc_assert, warning
from ffc.cpp import format
import numpy

def compute_integral_ir(itg_data,
                        form_data,
                        form_id,
                        parameters):
    "Compute intermediate represention of integral."

    info("Computing uflacs representation")

    # Initialise representation
    ir = initialize_integral_ir("uflacs", itg_data, form_data, form_id)

    import ffc.quadrature.quadraturerepresentation as qr
    import ffc.quadrature.quadratureutils as qu

    optimise_parameters = qr._parse_optimise_parameters(parameters)

    # Sort integrals and tabulate basis.
    sorted_integrals = qr._sort_integrals(itg_data.integrals, itg_data.metadata, form_data)
    integrals_dict, psi_tables, quad_weights = \
        qr._tabulate_basis(sorted_integrals, itg_data.domain_type, form_data)
    # TODO: Pass integrals_dict to uflacs, necessary for integral terms with different number of quadrature points
    ffc_assert(len(integrals_dict) == 1, "Assuming a single quadrature rule per integral domain for now.")

    # Save tables for quadrature weights and points
    ir["quadrature_weights"]  = quad_weights

    # Create dimensions of primary indices, needed to reset the argument 'A'
    # given to tabulate_tensor() by the assembler.
    ir["prim_idims"] = [qr.create_element(ufl_element).space_dimension()
                        for ufl_element in form_data.argument_elements]

    # FIXME: UFLACS: Understand how to use psi_tables
    (element_map, name_map, unique_tables) = \
        qu.create_psi_tables(psi_tables, optimise_parameters, vertex=None)

    from ffc.cpp import _generate_psi_name
    if 1:
        from pprint import pprint
        print '\n', '='*80, '\n'
        for ip, ip_tables in psi_tables.items():
            print '/'*25, 'ip:', ip
            for element, element_tables in ip_tables.items():
                print '_'*20, 'element:', str(element)
                element_number = element_map[ip][element]
                print 'element_number:', element_number
                component = ()
                for entity, entity_tables in element_tables.items():
                    print '-'*15, 'entity:', entity
                    for derivatives, values in entity_tables.items():
                        print '.'*10, 'derivatives:', derivatives
                        element_table_name = _generate_psi_name(element_number, entity, component, derivatives, vertex=None)
                        print 'element_table_name:', element_table_name
                        uvalues = unique_tables[element_table_name]
                        diff = numpy.transpose(uvalues) - values # Why is one transposed of the other? Which index is then ip and which is dofs?
                        pprint(values)
                        pprint(uvalues)
                        pprint(diff)
                        print numpy.sum(diff*diff)
        print '\n', '='*80, '\n', 'element_map:'
        pprint(element_map)
        print '\n', '='*80, '\n', 'name_map:'
        pprint(name_map)
        print '\n', '='*80, '\n', 'unique_tables:'
        pprint(unique_tables)
        print '\n', '='*80, '\n'

    # ------ End copy from quadraturerepresentation.py


    # TODO: Call upon ffc to build ir for element tables etc
    ffc_data = None


    # Delegate to flacs to build its intermediate representation and add to ir
    import uflacs.backends.ffc
    uir = uflacs.backends.ffc.compute_tabulate_tensor_ir(ffc_data, itg_data, form_data, parameters)
    ir.update(uir)

    return ir
