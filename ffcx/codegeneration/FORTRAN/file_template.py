# Code generation format strings for UFC (Unified Form-assembly Code)
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2018.

declaration_pre = """
! This code conforms with the UFC specification version {ufcx_version}
! and was automatically generated by FFCx version {ffcx_version}.
!
! This code was generated with the following options:
!
{options}

type :: ufcx_finite_element
  character, allocatable :: signature
  integer :: cell_shape
  integer :: element_type
  integer :: topological_dimension
  integer :: geometric_dimension
  integer :: space_dimension
  integer :: value_rank
  integer, dimension(:), allocatable :: value_shape
  integer :: value_size
  integer :: reference_value_rank
  integer, dimension(:), allocatable :: reference_value_shape
  integer :: reference_value_size
  integer :: degree
  integer :: block_size
  character, allocatable :: family
  integer :: basix_family
  integer :: basix_cell
  integer :: discontinuous
  integer :: lagrange_variant
  integer :: dpc_variant
  integer :: num_sub_elements
  integer :: sub_elements
  integer :: custom_element
end type

type :: ufcx_integral
   integer, dimension(:), allocatable :: enabled_coefficients
end type

type :: ufcx_form

end type

type :: ufcx_dofmap
  character, allocatable :: signature
  integer :: num_global_support_dofs
  integer :: num_element_support_dofs
  integer :: block_size
  integer, dimension(:), allocatable :: entity_dofs
  integer, dimension(:), allocatable :: entity_dof_offsets
  integer, dimension(:), allocatable :: entity_closure_dofs
  integer, dimension(:), allocatable :: entity_closure_dof_offsets
  integer :: num_sub_dofmaps
  type(ufcx_dofmap), pointer, dimension(:) :: sub_dofmaps
end type

"""

declaration_post = ""

implementation_pre = """
! This code conforms with the UFC specification version {ufcx_version}
! and was automatically generated by FFCx version {ffcx_version}.
!
! This code was generated with the following options:
!
{options}

include '{uflname}.f90h'
"""

implementation_post = """
"""
