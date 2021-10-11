/// This is UFC (Unified Form-assembly Code)
/// This code is released into the public domain.
///
/// The FEniCS Project (http://www.fenicsproject.org/) 2006-2019.
///
/// UFC defines the interface between code generated by FFCx and the
/// DOLFINx C++ library. Changes here must be reflected both in the FFCx
/// code generation and in the DOLFINx library calls.

#pragma once

#define UFC_VERSION_MAJOR 2021
#define UFC_VERSION_MINOR 1
#define UFC_VERSION_MAINTENANCE 0
#define UFC_VERSION_RELEASE 1

#if UFC_VERSION_RELEASE
#define UFC_VERSION                                                            \
  UFC_VERSION_MAJOR "." UFC_VERSION_MINOR "." UFC_VERSION_MAINTENANCE
#else
#define UFC_VERSION                                                            \
  UFC_VERSION_MAJOR "." UFC_VERSION_MINOR "." UFC_VERSION_MAINTENANCE ".dev0"
#endif

#include <stdbool.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C"
{

#if defined(__clang__)
#define restrict
#elif defined(__GNUC__) || defined(__GNUG__)
#define restrict __restrict__
#else
#define restrict
#endif // restrict
#endif // __cplusplus

  // <HEADER_DECL>

  typedef enum
  {
    interval = 10,
    triangle = 20,
    quadrilateral = 30,
    tetrahedron = 40,
    hexahedron = 50,
    vertex = 60,
    prism = 70,
    pyramid = 80
  } ufc_shape;

  typedef enum
  {
    cell = 0,
    exterior_facet = 1,
    interior_facet = 2
  } ufc_integral_type;

  typedef enum
  {
    ufc_basix_element = 0,
    ufc_mixed_element = 1,
    ufc_blocked_element = 2,
    ufc_quadrature_element = 3,
    ufc_custom_element = 4
  } ufc_element_type;

  /// Forward declarations
  typedef struct ufc_finite_element ufc_finite_element;
  typedef struct ufc_dofmap ufc_dofmap;

  // </HEADER_DECL>

  typedef struct ufc_finite_element
  {
    /// String identifying the finite element
    const char* signature;

    /// Return the cell shape
    ufc_shape cell_shape;

    /// Return the element type
    ufc_element_type element_type;

    /// Return the topological dimension of the cell shape
    int topological_dimension;

    /// Return the geometric dimension of the cell shape
    int geometric_dimension;

    /// Return the dimension of the finite element function space
    int space_dimension;

    /// Return the rank of the value space
    int value_rank;

    /// Return the dimension of the value space for axis i
    int* value_shape;

    /// Return the number of components of the value space
    int value_size;

    /// Return the rank of the reference value space
    int reference_value_rank;

    /// Return the dimension of the reference value space for axis i
    int* reference_value_shape;

    /// Return the number of components of the reference value space
    int reference_value_size;

    /// Return the maximum polynomial degree of the finite element
    /// function space
    int degree;

    /// Return the block size for a VectorElement. For a TensorElement,
    /// this is the product of the tensor's dimensions
    int block_size;

    /// Return the family of the finite element function space
    const char* family;

    /// Return the Basix identifier of the family of the finite element function space
    int basix_family;

    /// Return the Basix identifier of the cell shape
    int basix_cell;

    /// The Lagrange variant to be passed to Basix's create_element function
    int lagrange_variant;

    /// Indicates whether or not this is the discontinuous version of the element
    bool discontinuous;

    /// Return the number of sub elements (for a mixed element)
    int num_sub_elements;

    /// Get a finite element for sub element i (for a mixed
    /// element).
    ufc_finite_element** sub_elements;

  } ufc_finite_element;

  typedef struct ufc_dofmap
  {

    /// Return a string identifying the dofmap
    const char* signature;

    /// Number of dofs with global support (i.e. global constants)
    int num_global_support_dofs;

    /// Dimension of the local finite element function space for a cell
    /// (not including global support dofs)
    int num_element_support_dofs;

    /// Return the block size for a VectorElement or TensorElement
    int block_size;

    /// Number of dofs associated with each cell entity of dimension d
    int *num_entity_dofs;

    /// Tabulate the local-to-local mapping of dofs on entity (d, i)
    void (*tabulate_entity_dofs)(int* restrict dofs, int d, int i);

    /// Number of dofs associated with the closure of each cell entity of dimension d
    int *num_entity_closure_dofs;

    /// Tabulate the local-to-local mapping of dofs on the closure of entity (d, i)
    void (*tabulate_entity_closure_dofs)(int* restrict dofs, int d, int i);

    /// Return the number of sub dofmaps (for a mixed element)
    int num_sub_dofmaps;

    /// Get a dofmap for sub dofmap i (for a mixed element)
    ufc_dofmap** sub_dofmaps;

  } ufc_dofmap;

  /// Tabulate integral into tensor A with compiled quadrature rule
  ///
  /// @param[out] A
  /// @param[in] w Coefficients attached to the form to which the
  /// tabulated integral belongs.
  ///
  /// Dimensions: w[coefficient][restriction][dof].
  ///
  /// Restriction dimension
  /// applies to interior facet integrals, where coefficients restricted
  /// to both cells sharing the facet must be provided.
  /// @param[in] c Constants attached to the form to which the tabulated
  /// integral belongs. Dimensions: c[constant][dim].
  /// @param[in] coordinate_dofs Values of degrees of freedom of
  /// coordinate element. Defines the geometry of the cell. Dimensions:
  /// coordinate_dofs[restriction][num_dofs][3]. Restriction
  /// dimension applies to interior facet integrals, where cell
  /// geometries for both cells sharing the facet must be provided.
  /// @param[in] entity_local_index Local index of mesh entity on which
  /// to tabulate. This applies to facet integrals.
  /// @param[in] quadrature_permutation For facet integrals, numbers to
  /// indicate the permutation to be applied to each side of the facet
  /// to make the orientations of the faces matched up should be passed
  /// in. If an integer of value N is passed in, then:
  ///
  ///  - floor(N / 2) gives the number of rotations to apply to the
  ///  facet
  ///  - N % 2 gives the number of reflections to apply to the facet
  ///
  ///  For integrals not on facets, this argument has not effect and a
  ///  null pointer can be passed. For interior facets the array will
  ///  have size 2 (one permutation for each cell adjacent to the
  ///  facet). For exterior facets, this will have size 1.
  typedef void(ufc_tabulate_tensor_single)(
      float* restrict A, const float* restrict w,
      const float* restrict c, const double* restrict coordinate_dofs,
      const int* restrict entity_local_index,
      const uint8_t* restrict quadrature_permutation);

  /// Tabulate integral into tensor A with compiled
  /// quadrature rule and double precision
  ///
  /// @see ufc_tabulate_tensor_single
  ///
  typedef void(ufc_tabulate_tensor_double)(
      double* restrict A, const double* restrict w,
      const double* restrict c, const double* restrict coordinate_dofs,
      const int* restrict entity_local_index,
      const uint8_t* restrict quadrature_permutation);

  /// Tabulate integral into tensor A with compiled
  /// quadrature rule and complex single precision
  ///
  /// @see ufc_tabulate_tensor_single
  ///
  typedef void(ufc_tabulate_tensor_csingle)(
      float _Complex* restrict A, const float _Complex* restrict w,
      const float _Complex* restrict c, const double* restrict coordinate_dofs,
      const int* restrict entity_local_index,
      const uint8_t* restrict quadrature_permutation);

  /// Tabulate integral into tensor A with compiled
  /// quadrature rule and complex double precision
  ///
  /// @see ufc_tabulate_tensor_single
  ///
  typedef void(ufc_tabulate_tensor_cdouble)(
      double _Complex* restrict A, const double _Complex* restrict w,
      const double _Complex* restrict c, const double* restrict coordinate_dofs,
      const int* restrict entity_local_index,
      const uint8_t* restrict quadrature_permutation);

  typedef struct ufc_integral
  {
    const bool* enabled_coefficients;
    void* tabulate_tensor;
    bool needs_facet_permutations;
  } ufc_integral;

  typedef struct ufc_expression
  {

    /// Evaluate expression into tensor A with compiled evaluation
    /// points
    ///
    /// @param[out] A
    /// @param[in] w Coefficients attached to the expression.
    /// Dimensions: w[coefficient][dof].
    /// @param[in] c Constants attached to the expression. Dimensions:
    /// c[constant][dim].
    /// @param[in] coordinate_dofs Values of degrees of freedom of
    /// coordinate element. Defines the geometry of the cell.
    /// Dimensions: coordinate_dofs[num_dofs][3].
    void (*tabulate_expression)(double* restrict A,
                                const double* restrict w,
                                const double* restrict c,
                                const double* restrict coordinate_dofs);

    /// Positions of coefficients in original expression
    const int* original_coefficient_positions;

    /// Number of coefficients
    int num_coefficients;

    /// Number of evaluation points
    int num_points;

    /// Dimension of evaluation point, i.e. topological dimension of
    /// reference cell
    int topological_dimension;

    /// Indicates whether facet permutations are needed
    bool needs_facet_permutations;

    /// Coordinates of evaluations points. Dimensions:
    /// points[num_points][topological_dimension]
    const double* points;

    /// Return shape of expression. Dimension:
    /// value_shape[num_components]
    const int* value_shape;

    /// Number of components of return_shape
    int num_components;
  } ufc_expression;

  /// This class defines the interface for the assembly of the global
  /// tensor corresponding to a form with r + n arguments, that is, a
  /// mapping
  ///
  ///     a : V1 x V2 x ... Vr x W1 x W2 x ... x Wn -> R
  ///
  /// with arguments v1, v2, ..., vr, w1, w2, ..., wn. The rank r
  /// global tensor A is defined by
  ///
  ///     A = a(V1, V2, ..., Vr, w1, w2, ..., wn),
  ///
  /// where each argument Vj represents the application to the
  /// sequence of basis functions of Vj and w1, w2, ..., wn are given
  /// fixed functions (coefficients).
  typedef struct ufc_form
  {
    /// String identifying the form
    const char* signature;

    /// Rank of the global tensor (r)
    int rank;

    /// Number of coefficients (n)
    int num_coefficients;

    /// Number of constants
    int num_constants;

    /// Original coefficient position for each coefficient
    int* original_coefficient_position;

    /// Return list of names of coefficients
    const char** (*coefficient_name_map)(void);

    /// Return list of names of constants
    const char** (*constant_name_map)(void);

    /// Get a finite element for the i-th argument function, where 0 <=
    /// i < r + n.
    ///
    /// @param i Argument number if 0 <= i < r Coefficient number j = i
    /// - r if r + j <= i < r + n
    ufc_finite_element** finite_elements;

    /// Get a dofmap for the i-th argument function, where 0 <= i < r +
    /// n.
    ///
    /// @param i
    ///        Argument number if 0 <= i < r
    ///        Coefficient number j=i-r if r+j <= i < r+n
    ufc_dofmap** dofmaps;

    /// All ids for integrals
    int* (*integral_ids)(ufc_integral_type);

    /// Number of integrals
    int (*num_integrals)(ufc_integral_type);

    /// Get an integral on sub domain subdomain_id
    ufc_integral** (*integrals)(ufc_integral_type);

  } ufc_form;

  // FIXME: Formalise a UFC 'function space'
  typedef struct ufc_function_space
  {
    ufc_finite_element* finite_element;
    ufc_dofmap* dofmap;

    /// The family of the finite element for the geometry map
    const char* geometry_family;

    /// The degree of the finite element for the geometry map
    int geometry_degree;
  } ufc_function_space;

#ifdef __cplusplus
#undef restrict
}
#endif
