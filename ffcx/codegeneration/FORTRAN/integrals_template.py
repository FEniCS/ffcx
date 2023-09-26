# Code generation format strings for UFC (Unified Form-assembly Code)
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2018

declaration = """
type (ufcx_integral) :: {factory_name}
"""

factory = """
! Code for integral {factory_name}
contains
subroutine tabulate_tensor(A, w, c, coordinate_dofs, entity_local_index, quadrature_permutation)
 {scalar_type}, DIMENSION(0:*) :: A
 {scalar_type}, DIMENSION(0:*) :: w
 {scalar_type}, DIMENSION(0:*) :: c
 {geom_type}, DIMENSION(0:*) :: coordinate_dofs
 INTEGER, DIMENSION(0:*) :: entity_local_index
 CHARACTER, DIMENSION(0:*) :: quadrature_permutation
 {tabulate_tensor}
 contains
 function conditional(c, t, f)
   logical c
   double precision t, f, conditional
   if (c) then
     conditional = t
     return
   endif
   conditional = f
   return
 end function
 END subroutine

{factory_name}%enabled_coefficients = {enabled_coefficients}
{factory_name}%tabulate_tensor_{np_scalar_type} = tabulate_tensor_{factory_name}
{factory_name}%needs_facet_permutations = {needs_facet_permutations}
{factory_name}%coordinate_element = &{coordinate_element}
! End of code for integral {factory_name}
"""
