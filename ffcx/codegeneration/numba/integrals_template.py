# Code generation format strings for UFC (Unified Form-assembly Code)
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2018

factory = """
# Code for integral {factory_name}

c_signature = numba.types.void(
    numba.types.CPointer(numba.typeof(default_scalar_type())),
    numba.types.CPointer(numba.typeof(default_scalar_type())),
    numba.types.CPointer(numba.typeof(default_scalar_type())),
    numba.types.CPointer(numba.typeof(default_real_type())),
    numba.types.CPointer(numba.types.int32),
    numba.types.CPointer(numba.types.int32))


@numba.cfunc(c_signature, nopython=True)
def tabulate_tensor_{factory_name}(_A, _w, _c, _coordinate_dofs, _entity_local_index, _quadrature_permutation):
{tabulate_tensor}

class {factory_name}(object):
    enabled_coefficients = {enabled_coefficients}
    tabulate_tensor_{np_scalar_type} = tabulate_tensor_{factory_name}
    needs_facet_permutations = {needs_facet_permutations}
    coordinate_element = {coordinate_element}

# End of code for integral {factory_name}
"""
