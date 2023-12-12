/// This is UFCx
/// This software is released under the terms of the unlicense (see the file
/// UNLICENSE).
///
/// The FEniCS Project (http://www.fenicsproject.org/) 2006-2023.
///
/// UFCx defines the interface between code generated by FFCx and the
/// DOLFINx C++ library. Changes here must be reflected both in the FFCx
/// code generation and in the DOLFINx library calls.

#pragma once

#define UFCX_VERSION_MAJOR 0
#define UFCX_VERSION_MINOR 8
#define UFCX_VERSION_MAINTENANCE 0
#define UFCX_VERSION_RELEASE 0

#if UFCX_VERSION_RELEASE
#define UFCX_VERSION                                                           \
  UFCX_VERSION_MAJOR "." UFCX_VERSION_MINOR "." UFCX_VERSION_MAINTENANCE
#else
#define UFCX_VERSION                                                           \
  UFCX_VERSION_MAJOR "." UFCX_VERSION_MINOR "." UFCX_VERSION_MAINTENANCE ".de" \
                                                                         "v0"
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
  } ufcx_shape;

  typedef enum
  {
    cell = 0,
    exterior_facet = 1,
    interior_facet = 2
  } ufcx_integral_type;

  typedef enum
  {
    ufcx_basix_element = 0,
    ufcx_mixed_element = 1,
    ufcx_quadrature_element = 2,
    ufcx_basix_custom_element = 3,
    ufcx_real_element = 4,
  } ufcx_element_type;

  /// Forward declarations
  typedef struct ufcx_finite_element ufcx_finite_element;
  typedef struct ufcx_basix_custom_finite_element ufcx_basix_custom_finite_element;
  typedef struct ufcx_quadrature_rule ufcx_quadrature_rule;
  typedef struct ufcx_dofmap ufcx_dofmap;
  typedef struct ufcx_function_space ufcx_function_space;

  // </HEADER_DECL>

  typedef struct ufcx_finite_element
  {
    /// String identifying the finite element
    const char* signature;

    /// Cell shape
    ufcx_shape cell_shape;

    /// Element type
    ufcx_element_type element_type;

    /// Topological dimension of the cell
    int topological_dimension;

    /// Geometric dimension of the cell
    int geometric_dimension;

    /// Dimension of the finite element function space
    int space_dimension;

    /// Rank of the value space
    int value_rank;

    /// Dimension of the value space for axis i
    int* value_shape;

    /// Number of components of the value space
    int value_size;

    /// Rank of the reference value space
    int reference_value_rank;

    /// Dimension of the reference value space for axis i
    int* reference_value_shape;

    /// Number of components of the reference value space
    int reference_value_size;

    /// Maximum polynomial degree of the finite element function space
    int degree;

    /// Block size for a VectorElement. For a TensorElement, this is the
    /// product of the tensor's dimensions
    int block_size;

    /// Basix identifier of the family of the finite element function space
    int basix_family;

    /// Basix identifier of the cell shape
    int basix_cell;

    /// Indicates whether or not this is the discontinuous version of the
    /// element
    bool discontinuous;

    /// The Lagrange variant to be passed to Basix's create_element function
    int lagrange_variant;

    /// The DPC variant to be passed to Basix's create_element function
    int dpc_variant;

    /// Number of sub elements (for a mixed element)
    int num_sub_elements;

    /// Get a finite element for sub element i (for a mixed
    /// element).
    ufcx_finite_element** sub_elements;

    /// Pointer to data to recreate the element if it is a custom Basix element
    ufcx_basix_custom_finite_element* custom_element;

    /// Pointer to data to recreate the custom quadrature rule if the element has one
    ufcx_quadrature_rule* custom_quadrature;
  } ufcx_finite_element;

  typedef struct ufcx_basix_custom_finite_element
  {
    /// Basix identifier of the cell shape
    int cell_type;

    /// Dimension of the value space for axis i
    int value_shape_length;

    /// Dimension of the value space for axis i
    int* value_shape;

    /// The number of rows in the wcoeffs matrix
    int wcoeffs_rows;

    /// The number of columns in the wcoeffs matrix
    int wcoeffs_cols;

    /// The coefficients that define the polynomial set of the element in terms
    /// of the orthonormal polynomials on the cell
    double* wcoeffs;

    /// The number of interpolation points associated with each entity
    int* npts;

    /// The number of DOFs associated with each entity
    int* ndofs;

    // The coordinates of the interpolation points
    double* x;

    // The entries in the interpolation matrices
    double* M;

    /// The map type for the element
    int map_type;

    /// The Sobolev space for the element
    int sobolev_space;

    /// Indicates whether or not this is the discontinuous version of the
    /// element
    bool discontinuous;

    /// The highest degree full polynomial space contained in this element
    int embedded_subdegree;

    /// The number of derivatives needed when interpolating
    int interpolation_nderivs;

    /// The highest degree of a polynomial in the element
    int embedded_superdegree;

    /// The polyset type of the element
    int polyset_type;
  } ufcx_basix_custom_finite_element;

  typedef struct ufcx_quadrature_rule
  {
    /// Cell shape
    ufcx_shape cell_shape;

    /// The number of points
    int npts;

    /// The topological dimension of the cell
    int topological_dimension;

    /// The quadraute points
    double* points;

    /// The quadraute weights
    double* weights;
  } ufcx_quadrature_rule;

  typedef struct ufcx_dofmap
  {
    /// String identifying the dofmap
    const char* signature;

    /// Number of dofs with global support (i.e. global constants)
    int num_global_support_dofs;

    /// Dimension of the local finite element function space for a cell
    /// (not including global support dofs)
    int num_element_support_dofs;

    /// Return the block size for a VectorElement or TensorElement
    int block_size;

    /// Flattened list of dofs associated with each entity
    int* entity_dofs;

    /// Offset for dofs of each entity in entity_dofs
    int* entity_dof_offsets;

    /// Flattened list of closure dofs associated with each entity
    int* entity_closure_dofs;

    /// Offset for closure dofs of each entity in entity_closure_dofs
    int* entity_closure_dof_offsets;

    /// Number of sub dofmaps (for a mixed element)
    int num_sub_dofmaps;

    /// Get a dofmap for sub dofmap i (for a mixed element)
    ufcx_dofmap** sub_dofmaps;

  } ufcx_dofmap;

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
  /// For integrals not on interior facets, this argument has no effect and a
  /// null pointer can be passed. For interior facets the array will have size 2
  /// (one permutation for each cell adjacent to the facet).
  typedef void(ufcx_tabulate_tensor_float32)(
      float* restrict A, const float* restrict w, const float* restrict c,
      const float* restrict coordinate_dofs,
      const int* restrict entity_local_index,
      const uint8_t* restrict quadrature_permutation);

  /// Tabulate integral into tensor A with compiled
  /// quadrature rule and double precision
  ///
  /// @see ufcx_tabulate_tensor_single
  typedef void(ufcx_tabulate_tensor_float64)(
      double* restrict A, const double* restrict w, const double* restrict c,
      const double* restrict coordinate_dofs,
      const int* restrict entity_local_index,
      const uint8_t* restrict quadrature_permutation);

  /// Tabulate integral into tensor A with compiled
  /// quadrature rule and extended double precision
  ///
  /// @see ufcx_tabulate_tensor_single
  typedef void(ufcx_tabulate_tensor_longdouble)(
      long double* restrict A, const long double* restrict w,
      const long double* restrict c,
      const long double* restrict coordinate_dofs,
      const int* restrict entity_local_index,
      const uint8_t* restrict quadrature_permutation);

  /// Tabulate integral into tensor A with compiled
  /// quadrature rule and complex single precision
  ///
  /// @see ufcx_tabulate_tensor_single
  typedef void(ufcx_tabulate_tensor_complex64)(
      float _Complex* restrict A, const float _Complex* restrict w,
      const float _Complex* restrict c, const float* restrict coordinate_dofs,
      const int* restrict entity_local_index,
      const uint8_t* restrict quadrature_permutation);

  /// Tabulate integral into tensor A with compiled
  /// quadrature rule and complex double precision
  ///
  /// @see ufcx_tabulate_tensor_single
  typedef void(ufcx_tabulate_tensor_complex128)(
      double _Complex* restrict A, const double _Complex* restrict w,
      const double _Complex* restrict c, const double* restrict coordinate_dofs,
      const int* restrict entity_local_index,
      const uint8_t* restrict quadrature_permutation);

  typedef struct ufcx_integral
  {
    const bool* enabled_coefficients;
    ufcx_tabulate_tensor_float32* tabulate_tensor_float32;
    ufcx_tabulate_tensor_float64* tabulate_tensor_float64;
    ufcx_tabulate_tensor_longdouble* tabulate_tensor_longdouble;
    ufcx_tabulate_tensor_complex64* tabulate_tensor_complex64;
    ufcx_tabulate_tensor_complex128* tabulate_tensor_complex128;
    bool needs_facet_permutations;

    /// Get the coordinate element associated with the geometry of the mesh.
    ufcx_finite_element* coordinate_element;
  } ufcx_integral;

  typedef struct ufcx_expression
  {
    /// Evaluate expression into tensor A with compiled evaluation points
    ///
    /// @param[out] A
    ///         Dimensions: A[num_points][num_components][num_argument_dofs]
    ///
    /// @see ufcx_tabulate_tensor
    ///
    ufcx_tabulate_tensor_float32* tabulate_tensor_float32;
    ufcx_tabulate_tensor_float64* tabulate_tensor_float64;
    ufcx_tabulate_tensor_longdouble* tabulate_tensor_longdouble;
    ufcx_tabulate_tensor_complex64* tabulate_tensor_complex64;
    ufcx_tabulate_tensor_complex128* tabulate_tensor_complex128;

    /// Number of coefficients
    int num_coefficients;

    /// Number of constants
    int num_constants;

    /// Original coefficient position for each coefficient
    const int* original_coefficient_positions;

    /// List of names of coefficients
    const char** coefficient_names;

    /// List of names of constants
    const char** constant_names;

    /// Number of evaluation points
    int num_points;

    /// Dimension of evaluation point, i.e. topological dimension of
    /// reference cell
    int topological_dimension;

    /// Coordinates of evaluations points. Dimensions:
    /// points[num_points][topological_dimension]
    const double* points;

    /// Shape of expression. Dimension: value_shape[num_components]
    const int* value_shape;

    /// Number of components of return_shape
    int num_components;

    /// Rank, i.e. number of arguments
    int rank;

    /// Function spaces for all functions in the Expression.
    ///
    /// Function spaces for coefficients are followed by
    /// Arguments function spaces.
    /// Dimensions: function_spaces[num_coefficients + rank]
    ufcx_function_space** function_spaces;

  } ufcx_expression;

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
  typedef struct ufcx_form
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

    /// List of names of coefficients
    const char** coefficient_name_map;

    /// List of names of constants
    const char** constant_name_map;

    /// Get a finite element for the i-th argument function, where 0 <=
    /// i < r + n.
    ///
    /// @param i Argument number if 0 <= i < r Coefficient number j = i
    /// - r if r + j <= i < r + n
    ufcx_finite_element** finite_elements;

    /// Get a dofmap for the i-th argument function, where 0 <= i < r +
    /// n.
    ///
    /// @param i
    ///        Argument number if 0 <= i < r
    ///        Coefficient number j=i-r if r+j <= i < r+n
    ufcx_dofmap** dofmaps;

    /// List of cell, interior facet and exterior facet integrals
    ufcx_integral** form_integrals;

    /// IDs for each integral in form_integrals list
    int* form_integral_ids;

    /// Offsets for cell, interior facet and exterior facet integrals in
    /// form_integrals list
    int* form_integral_offsets;

  } ufcx_form;

  // FIXME: Formalise a UFCX 'function space'
  typedef struct ufcx_function_space
  {
    ufcx_finite_element* finite_element;
    ufcx_dofmap* dofmap;

    /// The family of the finite element for the geometry map
    const char* geometry_family;

    /// The degree of the finite element for the geometry map
    int geometry_degree;

    /// The Basix cell of the finite element for the geometry map
    int geometry_basix_cell;

    /// The Basix variant of the finite element for the geometry map
    int geometry_basix_variant;
  } ufcx_function_space;

#ifdef __cplusplus
#undef restrict
}
#endif
