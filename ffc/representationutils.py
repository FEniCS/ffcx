"""This module contains utility functions for some code shared between
quadrature and tensor representation."""

from ffc.fiatinterface import create_element

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

    # Extract reference and physical value dimensions
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

    # Figure out which sub-element number we are in:
    tot = physical_value_dims[0]
    for (sub_element_number, i) in enumerate(physical_value_dims):
        if component < tot:
            break
        else:
            tot += i

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
