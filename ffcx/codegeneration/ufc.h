/// This is UFC (Unified Form-assembly Code)
/// This code is released into the public domain.
///
/// The FEniCS Project (http://www.fenicsproject.org/) 2006-2019.
///
/// UFC defines the interface between code generated by FFCX and the
/// DOLFIN C++ library. Changes here must be reflected both in the FFCX
/// code generation and in the DOLFIN library calls.

#pragma once

#define UFC_VERSION_MAJOR 2018
#define UFC_VERSION_MINOR 1
#define UFC_VERSION_MAINTENANCE 0
#define UFC_VERSION_RELEASE 0

#if UFC_VERSION_RELEASE
#define UFC_VERSION                                                            \
  UFC_VERSION_MAJOR "." UFC_VERSION_MINOR "." UFC_VERSION_MAINTENANCE
#else
#define UFC_VERSION                                                            \
  UFC_VERSION_MAJOR "." UFC_VERSION_MINOR "." UFC_VERSION_MAINTENANCE ".dev0"
#endif

#include <stdbool.h>
#include <stdint.h>
#include <ufc_geometry.h>

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
  } ufc_shape;

  typedef enum
  {
    PointEval = 1,
    ComponentPointEval = 2,
    PointNormalDeriv = 3,
    IntegralMoment = 4,
    FrobeniusIntegralMoment = 5,
    PointEdgeTangent = 6,
    PointFaceTangent = 7,
    PointScaledNormalEval = 8,
    PointDeriv = 9,
    IntegralMomentOfNormalDerivative = 10,
    PointNormalEval = 11,
    PointwiseInnerProductEval = 12,
  } ufc_doftype;

  /// Forward declarations
  typedef struct ufc_coordinate_mapping ufc_coordinate_mapping;
  typedef struct ufc_finite_element ufc_finite_element;
  typedef struct ufc_dofmap ufc_dofmap;

  // </HEADER_DECL>

  typedef struct ufc_finite_element
  {
    /// String identifying the finite element
    const char* signature;

    /// Return the cell shape
    ufc_shape cell_shape;

    /// Return the topological dimension of the cell shape
    int topological_dimension;

    /// Return the geometric dimension of the cell shape
    int geometric_dimension;

    /// Return the dimension of the finite element function space
    int space_dimension;

    /// Return the rank of the value space
    int value_rank;

    /// Return the dimension of the value space for axis i
    int (*value_dimension)(int i);

    /// Return the number of components of the value space
    int value_size;

    /// Return the rank of the reference value space
    int reference_value_rank;

    /// Return the dimension of the reference value space for axis i
    int (*reference_value_dimension)(int i);

    /// Return the number of components of the reference value space
    int reference_value_size;

    /// Return the maximum polynomial degree of the finite element
    /// function space
    int degree;

    /// Return the block size for a VectorElement.
    /// For a TensorElement, this is the product of the tensor's dimensions
    int block_size;

    /// Return the family of the finite element function space
    const char* family;

    int (*evaluate_reference_basis)(double* restrict reference_values,
                                    int num_points, const double* restrict X);

    int (*evaluate_reference_basis_derivatives)(
        double* restrict reference_values, int order, int num_points,
        const double* restrict X);

    /// Map the reference values on to the cell
    /// @param[in] cell_permutations An integer that says how each entity of the
    ///         cell of dimension < tdim has been permuted relative to a
    ///         low-to-high ordering of the cell.

    int (*transform_reference_basis_derivatives)(
        double* restrict values, int order, int num_points,
        const double* restrict reference_values, const double* restrict X,
        const double* restrict J, const double* restrict detJ,
        const double* restrict K, const uint32_t cell_permutation);

    /// Map values of field from physical to reference space which has
    /// been evaluated at points given by
    /// tabulate_reference_dof_coordinates.
    int (*transform_values)(ufc_scalar_t* restrict reference_values,
                            const ufc_scalar_t* restrict physical_values,
                            const double* restrict coordinate_dofs,
                            const ufc_coordinate_mapping* restrict cm);

    // FIXME: change to 'const double* reference_dof_coordinates()'
    /// Tabulate the coordinates of all dofs on a reference cell
    int (*tabulate_reference_dof_coordinates)(
        double* restrict reference_dof_coordinates);

    /// Return the number of sub elements (for a mixed element)
    int num_sub_elements;

    /// Create a new finite element for sub element i (for a mixed
    /// element). Memory for the new object is obtained with malloc(),
    /// and the caller is reponsible for freeing it by calling free().
    ufc_finite_element* (*create_sub_element)(int i);

    /// Create a new class instance. Memory for the new object is
    /// obtained with malloc(), and the caller is reponsible for
    /// freeing it by calling free().
    ufc_finite_element* (*create)(void);
  } ufc_finite_element;

  typedef struct ufc_dofmap
  {

    /// Return a string identifying the dofmap
    const char* signature;

    /// The base DoF permutations
    const int* base_permutations;

    /// The size of the base DoF permutations array
    int size_base_permutations;

    /// Number of dofs with global support (i.e. global constants)
    int num_global_support_dofs;

    /// Dimension of the local finite element function space
    /// for a cell (not including global support dofs)
    int num_element_support_dofs;

    /// Return the block size for a VectorElement
    int submap_block_size;

    /// Number of dofs associated with each cell entity of
    /// dimension d
    int num_entity_dofs[4];

    /// Tabulate the local-to-local mapping of dofs on entity (d, i)
    void (*tabulate_entity_dofs)(int* restrict dofs, int d, int i);

    /// Return the number of sub dofmaps (for a mixed element)
    int num_sub_dofmaps;

    /// Create a new dofmap for sub dofmap i (for a mixed
    /// element). Memory for the new object is obtained with malloc(),
    /// and the caller is reponsible for freeing it by calling free().
    ufc_dofmap* (*create_sub_dofmap)(int i);

    /// Create a new class instance. Memory for the new object is
    /// obtained with malloc(), and the caller is reponsible for
    /// freeing it by calling free().
    ufc_dofmap* (*create)(void);
  } ufc_dofmap;

  /// A representation of a coordinate mapping parameterized by a local
  /// finite element basis on each cell
  typedef struct ufc_coordinate_mapping
  {

    /// Return coordinate_mapping signature string
    const char* signature;

    /// Create object of the same type. Memory for the new object is
    /// obtained with malloc(), and the caller is reponsible for
    /// freeing it by calling free().
    ufc_coordinate_mapping* (*create)(void);

    /// Return geometric dimension of the coordinate_mapping
    int geometric_dimension;

    /// Return topological dimension of the coordinate_mapping
    int topological_dimension;

    /// Boolean flag for affine 
    int is_affine;

    /// Return cell shape of the coordinate_mapping
    ufc_shape cell_shape;

    /// Create dofmap for the underlying scalar element. Memory for
    /// the new object is obtained with malloc(), and the caller is
    /// reponsible for freeing it by calling free().
    ufc_dofmap* (*create_scalar_dofmap)(void);

    /// Evaluate basis (and derivatives) of the associated element
    int (*evaluate_basis_derivatives)(
        double* restrict reference_values, int order, int num_points,
        const double* restrict X);

  } ufc_coordinate_mapping;

  /// Tabulate integral into tensor A with compiled quadrature rule
  ///
  /// @param[out] A
  /// @param[in] w Coefficients attached to the form to which the
  ///         tabulated integral belongs.
  ///
  ///         Dimensions: w[coefficient][restriction][dof].
  ///
  ///         Restriction dimension
  ///         applies to interior facet integrals, where coefficients
  ///         restricted to both cells sharing the facet must be
  ///         provided.
  /// @param[in] c Constants attached to the form to which the tabulated
  ///         integral belongs. Dimensions: c[constant][dim].
  /// @param[in] coordinate_dofs Values of degrees of freedom of
  ///         coordinate element. Defines the geometry of the cell.
  ///         Dimensions: coordinate_dofs[restriction][num_dofs][gdim].
  ///         Restriction dimension applies to interior facet integrals,
  ///         where cell geometries for both cells sharing the facet
  ///         must be provided.
  /// @param[in] entity_local_index Local index of mesh entity on which
  ///         to tabulate. This applies to facet integrals.
  /// @param[in] quadrature_permutation For facet integrals, numbers to
  ///         indicate the permutation to be applied to each side of the
  ///         facet to make the orientations of the faces matched up
  ///         should be passed in. If an integer of value N is passed
  ///         in, then:
  ///
  ///          - floor(N / 2) gives the number of rotations to apply to the
  ///          facet
  ///          - N % 2 gives the number of reflections to apply to the facet
  ///
  ///         For integrals not on facets, this argument has not effect
  ///         and a null pointer can be passed. For
  ///         interior facets the array will have size 2 (one permutation
  ///         for each cell adjacent to the facet). For exterior facets,
  ///         this will have size 1.
  /// @param[in] cell_permutations An integer that says how each entity of the
  ///         cell of dimension < tdim has been permuted relative to a
  ///         low-to-high ordering of the cell. This bits of this integer
  ///         represent (from least to most significant bits):
  ///
  ///          - Faces (3 bits each). Reflections are least significant bit,
  ///          then next two bits give number of rotations.
  ///          - Edges (1 bit each). The bit is 1 if the edge is reflected.
  ///
  typedef void(ufc_tabulate_tensor)(
      ufc_scalar_t* restrict A, const ufc_scalar_t* restrict w,
      const ufc_scalar_t* restrict c, const double* restrict coordinate_dofs,
      const int* restrict entity_local_index,
      const uint8_t* restrict quadrature_permutation,
      const uint32_t cell_permutation);

  /// Tabulate integral into tensor A with runtime quadrature rule
  ///
  /// @see ufc_tabulate_tensor
  ///
  typedef void(ufc_tabulate_tensor_custom)(
      ufc_scalar_t* restrict A, const ufc_scalar_t* restrict w,
      const ufc_scalar_t* restrict c, const double* restrict coordinate_dofs,
      int num_quadrature_points, const double* restrict quadrature_points,
      const double* restrict quadrature_weights,
      const double* restrict facet_normals);

  typedef struct ufc_integral
  {
    const bool* enabled_coefficients;
    ufc_tabulate_tensor* tabulate_tensor;
  } ufc_integral;

  typedef struct ufc_custom_integral
  {
    const bool* enabled_coefficients;
    ufc_tabulate_tensor_custom* tabulate_tensor;
  } ufc_custom_integral;

  typedef struct ufc_expression
  {

    /// Evaluate expression into tensor A with compiled evaluation points
    ///
    /// @param[out] A
    /// @param[in] w
    ///         Coefficients attached to the expression.
    ///         Dimensions: w[coefficient][dof].
    /// @param[in] c
    ///         Constants attached to the expression.
    ///         Dimensions: c[constant][dim].
    /// @param[in] coordinate_dofs
    ///         Values of degrees of freedom of coordinate element.
    ///         Defines the geometry of the cell.
    ///         Dimensions: coordinate_dofs[num_dofs][gdim].
    ///
    void (*tabulate_expression)(ufc_scalar_t* restrict A,
                                const ufc_scalar_t* restrict w,
                                const ufc_scalar_t* restrict c,
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

    /// Coordinates of evaluations points.
    /// Dimensions: points[num_points][topological_dimension].
    const double* points;

    /// Return shape of expression
    /// Dimension: value_shape[num_components].
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

    /// Return original coefficient position for each coefficient
    ///
    /// @param i
    ///        Coefficient number, 0 <= i < n
    ///
    int (*original_coefficient_position)(int i);

    /// Return list of names of coefficients
    const char** (*coefficient_name_map)(void);

    /// Return list of names of constants
    const char** (*constant_name_map)(void);

    /// Create a new coordinate mapping. Memory for the new object is
    /// obtained with malloc(), and the caller is reponsible for
    /// freeing it by calling free().
    ufc_coordinate_mapping* (*create_coordinate_mapping)(void);

    /// Create a new finite element for the i-th argument function,
    /// where 0 <= i < r+n. Memory for the new object is obtained with
    /// malloc(), and the caller is reponsible for freeing it by
    /// calling free().
    ///
    /// @param i
    ///        Argument number if 0 <= i < r
    ///        Coefficient number j=i-r if r+j <= i < r+n
    ///
    ufc_finite_element* (*create_finite_element)(int i);

    /// Create a new dofmap for the i-th argument function, where 0 <=
    /// i < r+n.  Memory for the new object is obtained with malloc(),
    /// and the caller is reponsible for freeing it by calling free().
    ///
    /// @param i
    ///        Argument number if 0 <= i < r
    ///        Coefficient number j=i-r if r+j <= i < r+n
    ///
    ufc_dofmap* (*create_dofmap)(int i);

    /// All ids for cell integrals
    void (*get_cell_integral_ids)(int* ids);

    /// All ids for exterior facet integrals
    void (*get_exterior_facet_integral_ids)(int* ids);

    /// All ids for interior facet integrals
    void (*get_interior_facet_integral_ids)(int* ids);

    /// All ids for vertex integrals
    void (*get_vertex_integral_ids)(int* ids);

    /// All ids for custom integrals
    void (*get_custom_integral_ids)(int* ids);

    /// Number of cell integrals
    int num_cell_integrals;

    /// Number of exterior facet integrals
    int num_exterior_facet_integrals;

    /// Number of interior facet integrals
    int num_interior_facet_integrals;

    /// Number of vertex integrals
    int num_vertex_integrals;

    /// Number of custom integrals
    int num_custom_integrals;

    /// Create a new cell integral on sub domain subdomain_id. Memory
    /// for the new object is obtained with malloc(), and the caller
    /// is reponsible for freeing it by calling free().
    ufc_integral* (*create_cell_integral)(int subdomain_id);

    /// Create a new exterior facet integral on sub domain
    /// subdomain_id. Memory for the new object is obtained with
    /// malloc(), and the caller is reponsible for freeing it by
    /// calling free().
    ufc_integral* (*create_exterior_facet_integral)(int subdomain_id);

    /// Create a new interior facet integral on sub domain
    /// subdomain_id. Memory for the new object is obtained with
    /// malloc(), and the caller is reponsible for freeing it by
    /// calling free().
    ufc_integral* (*create_interior_facet_integral)(int subdomain_id);

    /// Create a new vertex integral on sub domain
    /// subdomain_id. Memory for the new object is obtained with
    /// malloc(), and the caller is reponsible for freeing it by
    /// calling free().
    ufc_integral* (*create_vertex_integral)(int subdomain_id);

    /// Create a new custom integral on sub domain
    /// subdomain_id. Memory for the new object is obtained with
    /// malloc(), and the caller is reponsible for freeing it by
    /// calling free().
    ufc_custom_integral* (*create_custom_integral)(int subdomain_id);

  } ufc_form;

  // FIXME: Formalise a UFC 'function space'.
  typedef struct ufc_function_space
  {
    // Pointer to factory function that creates a new
    // ufc_finite_element. Memory for the new object is obtained with
    // malloc(), and the caller is reponsible for freeing it by
    // calling free().
    ufc_finite_element* (*create_element)(void);

    // Pointer to factory function that creates a new
    // ufc_dofmap. Memory for the new object is obtained with
    // malloc(), and the caller is reponsible for freeing it by
    // calling free().
    ufc_dofmap* (*create_dofmap)(void);

    // Pointer to factory function that creates a new
    // ufc_coordinate_mapping. Memory for the new object is obtained
    // with malloc(), and the caller is reponsible for freeing it by
    // calling free().
    ufc_coordinate_mapping* (*create_coordinate_mapping)(void);
  } ufc_function_space;

#ifdef __cplusplus
#undef restrict
}
#endif
