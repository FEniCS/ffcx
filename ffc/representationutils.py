"""This module contains utility functions for some code shared between
quadrature and tensor representation."""

# Copyright (C) 2012-2013 Marie Rognes
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
# Modified by Martin Alnaes, 2013
#
# First added:  2013-01-08
# Last changed: 2013-02-10

from ffc.fiatinterface import create_element
from ffc.fiatinterface import cellname_to_num_entities
from ffc.cpp import format

def transform_component(component, offset, ufl_element):
    """
    This function accounts for the fact that if the geometrical and
    topological dimension does not match, then for native vector
    elements, in particular the Piola-mapped ones, the physical value
    dimensions and the reference value dimensions are not the
    same. This has certain consequences for mixed elements, aka 'fun
    with offsets'.
    """
    # This code is used for tensor/monomialtransformation.py and
    # quadrature/quadraturetransformerbase.py.

    gdim = ufl_element.cell().geometric_dimension()
    tdim = ufl_element.cell().topological_dimension()

    # Do nothing if we are not in a special case: The special cases
    # occur if we have piola mapped elements (for which value_shape !=
    # ()), and if gdim != tdim)
    if gdim == tdim:
        return component, offset
    all_mappings =  create_element(ufl_element).mapping()
    special_case = (any(['piola' in m for m in all_mappings])
                    and ufl_element.num_sub_elements() > 1)
    if not special_case:
        return component, offset

    # Extract lists of reference and physical value dimensions by
    # sub-element
    reference_value_dims = []
    physical_value_dims = []
    for sub_element in ufl_element.sub_elements():
        assert (len(sub_element.value_shape()) < 2), \
            "Vector-valued assumption failed"
        if sub_element.value_shape() == ():
            reference_value_dims += [1]
            physical_value_dims += [1]
        else:
            reference_value_dims += [sub_element.value_shape()[0]
                                     - (gdim - tdim)]
            physical_value_dims += [sub_element.value_shape()[0]]

    # Figure out which sub-element number 'component' is in,
    # 'sub_element_number' contains the result
    tot = physical_value_dims[0]
    for sub_element_number in range(len(physical_value_dims)):
        if component < tot:
            break
        else:
            tot += physical_value_dims[sub_element_number+1]

    # Compute the new reference offset:
    reference_offset = sum(reference_value_dims[:sub_element_number])
    physical_offset = sum(physical_value_dims[:sub_element_number])
    shift = physical_offset - reference_offset

    # Compute the component relative to the reference frame
    reference_component = component - shift

    return reference_component, reference_offset

def needs_oriented_jacobian(form_data):
    # Check whether this form needs an oriented jacobian (only forms
    # involgin contravariant piola mappings seem to need it)
    for ufl_element in form_data.unique_elements:
        element = create_element(ufl_element)
        if "contravariant piola" in element.mapping():
            return True
    return False

def initialize_integral_ir(representation, itg_data, form_data, form_id):
    """Initialize a representation dict with common information that is
    expected independently of which representation is chosen."""
    entitytype = { "cell": "cell",
                   "exterior_facet": "facet",
                   "interior_facet": "facet",
                   "point": "vertex",
                   }[itg_data.domain_type]
    return { "representation":       representation,
             "domain_type":          itg_data.domain_type,
             "domain_id":            itg_data.domain_id,
             "form_id":              form_id,
             "rank":                 form_data.rank,
             "geometric_dimension":  form_data.geometric_dimension,
             "topological_dimension":form_data.topological_dimension,
             "entitytype":           entitytype,
             "num_facets":           cellname_to_num_entities[form_data.cell.cellname()][-2],
             "num_vertices":         cellname_to_num_entities[form_data.cell.cellname()][0],
             "needs_oriented":       needs_oriented_jacobian(form_data),
           }

def initialize_integral_code(ir, prefix, parameters):
    "Representation independent default initialization of code dict for integral from intermediate representation."
    code = {}
    code["restrict"] = parameters["restrict_keyword"]
    code["classname"] = format["classname " + ir["domain_type"] + "_integral"](prefix, ir["form_id"], ir["domain_id"])
    code["members"] = ""
    code["constructor"] = format["do nothing"]
    code["constructor_arguments"] = ""
    code["initializer_list"] = ""
    code["destructor"] = format["do nothing"]
    #code["additional_includes_set"] = set() #ir["additional_includes_set"]
    return code
