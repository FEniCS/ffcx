# Copyright (C) 2025 Chris Richardson and Paul T. Kühner
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Template for integral output."""

factory = """
# Code for integral {factory_name}

def tabulate_tensor_{factory_name}(_A, _w, _c, _coordinate_dofs,
                                   _entity_local_index, _quadrature_permutation, custom_data):
{tabulate_tensor}

class {factory_name}(object):
    enabled_coefficients = {enabled_coefficients}
    tabulate_tensor = tabulate_tensor_{factory_name}
    needs_facet_permutations = {needs_facet_permutations}
    coordinate_element_hash = {coordinate_element_hash}
    domain = {domain}

# End of code for integral {factory_name}
"""
