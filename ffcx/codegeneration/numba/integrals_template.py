# Code generation format strings for UFC (Unified Form-assembly Code)
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2018

factory = """
# Code for integral {factory_name}
import numba

def tabulate_tensor_{factory_name}(_A, _w, _c, _coordinate_dofs, _entity_local_index, _quadrature_permutation):
{tabulate_tensor}

class {factory_name}(object):
    enabled_coefficients = {enabled_coefficients}
    tabulate_tensor = tabulate_tensor_{factory_name}
    needs_facet_permutations = {needs_facet_permutations}
    coordinate_element = {coordinate_element}

# End of code for integral {factory_name}
"""
