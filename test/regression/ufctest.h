// Copyright (C) 2010 Anders Logg.
// Licensed under the GNU GPL version 3 or any later version.
//
// First added:  2010-01-24
// Last changed: 2010-01-28
//
// Functions for calling generated UFC functions with "random" (but
// fixed) data and print the output to screen. Useful for running
// regression tests.

#include <sstream>
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <ufc.h>

typedef unsigned int uint;

// How many derivatives to test
const uint max_derivative = 2;

// Precision in output of floats
const uint precision = 10;

// Global counter for results
uint counter = 0;

// Function for printing a single value
template <class value_type>
void print_value(value_type value)
{
  std::cout.precision(precision);
  if (std::abs(static_cast<double>(value)) < std::pow(0.1, precision))
    std::cout << "0";
  else
    std::cout << value;
}

// Function for printing scalar result
template <class value_type>
void print_scalar(std::string name, value_type value, int i=-1, int j=-1)
{
  std::stringstream s;
  s << counter++ << "_";
  s << name;
  if (i >= 0) s << "_" << i;
  if (j >= 0) s << "_" << j;
  std::cout << s.str() << " = ";
  print_value(value);
  std::cout << std::endl;
}

// Function for printing array result
template <class value_type>
void print_array(std::string name, unsigned int n, value_type* values, int i=-1, int j=-1)
{
  std::stringstream s;
  s << counter++ << "_";
  s << name;
  if (i >= 0) s << "_" << i;
  if (j >= 0) s << "_" << j;
  std::cout << s.str() << " =";
  for (uint i = 0; i < n; i++)
  {
    std::cout << " ";
    print_value(values[i]);
  }
  std::cout << std::endl;
}

// Class for creating a "random" ufc::mesh object
class test_mesh : public ufc::mesh
{
public:

  test_mesh(ufc::shape cell_shape)
  {
    // Store dimensions
    switch (cell_shape)
    {
    case ufc::interval:
      topological_dimension = 1;
      geometric_dimension = 1;
      break;
    case ufc::triangle:
      topological_dimension = 2;
      geometric_dimension = 2;
      break;
    case ufc::tetrahedron:
      topological_dimension = 3;
      geometric_dimension = 3;
      break;
    default:
      throw std::runtime_error("Unhandled cell shape.");
    }

    // Set some random sizes
    num_entities = new uint[4];
    num_entities[0] = 10001;
    num_entities[1] = 10002;
    num_entities[2] = 10003;
    num_entities[3] = 10004;
  }

  ~test_mesh()
  {
    delete [] num_entities;
  }

};

// Class for creating "random" ufc::cell objects
class test_cell : public ufc::cell
{
public:

  test_cell(ufc::shape cell_shape, int offset=0)
  {
    // Store cell shape
    this->cell_shape = cell_shape;

    // Store dimensions
    switch (cell_shape)
    {
    case ufc::interval:
      topological_dimension = 1;
      geometric_dimension = 1;
      break;
    case ufc::triangle:
      topological_dimension = 2;
      geometric_dimension = 2;
      break;
    case ufc::tetrahedron:
      topological_dimension = 3;
      geometric_dimension = 3;
      break;
    default:
      throw std::runtime_error("Unhandled cell shape.");
    }

    // Generate some "random" entity indices
    entity_indices = new uint * [4];
    for (uint i = 0; i < 4; i++)
    {
      entity_indices[i] = new uint[6];
      for (uint j = 0; j < 6; j++)
        entity_indices[i][j] = i*j + offset;
    }

    // Generate some "random" coordinates
    double** x = new double * [4];
    for (uint i = 0; i < 4; i++)
      x[i] = new double[3];
    x[0][0] = 0.90; x[0][1] = 0.34; x[0][2] = 0.45;
    x[1][0] = 0.56; x[1][1] = 0.76; x[1][2] = 0.83;
    x[2][0] = 0.98; x[2][1] = 0.78; x[2][2] = 0.19;
    x[3][0] = 0.12; x[3][1] = 0.56; x[3][2] = 0.66;
    coordinates = x;
  }

  ~test_cell()
  {
    for (uint i = 0; i < 4; i++)
    {
      delete [] entity_indices[i];
      delete [] coordinates[i];
    }
    delete [] entity_indices;
    delete [] coordinates;
  }

};

// Class for creating a "random" ufc::function object
class test_function : public ufc::function
{
public:

  test_function(uint value_size) : value_size(value_size) {}

  void evaluate(double* values, const double* coordinates, const ufc::cell& c) const
  {
    for (uint i = 0; i < value_size; i++)
    {
      values[i] = 1.0;
      for (uint j = 0; j < c.geometric_dimension; j++)
        values[i] *= static_cast<double>(i + 1)*coordinates[j];
    }
  }

private:

  uint value_size;

};

// Function for testing ufc::element objects
void test_finite_element(ufc::finite_element& element)
{
  std::cout << std::endl;
  std::cout << "Testing finite_element" << std::endl;
  std::cout << "----------------------" << std::endl;

  // Prepare arguments
  test_cell c(element.cell_shape());
  uint value_size = 1;
  for (uint i = 0; i < element.value_rank(); i++)
    value_size *= element.value_dimension(i);
  uint derivative_size = 1;
  for (uint i = 0; i < max_derivative; i++)
    derivative_size *= c.geometric_dimension;
  double* values = new double[element.space_dimension()*value_size*derivative_size];
  for (uint i = 0; i < element.space_dimension()*value_size*derivative_size; i++)
  {
    values[i] = 0.0;
  }
  double* dof_values = new double[element.space_dimension()];
  for (uint i = 0; i < element.space_dimension(); i++)
  {
    dof_values[i] = 0.0;
  }
  double* vertex_values = new double[(c.topological_dimension + 1)*value_size];
  for (uint i = 0; i < (c.topological_dimension + 1)*value_size; i++)
  {
    vertex_values[i] = 0.0;
  }
  double* coordinates = new double[c.geometric_dimension];
  for (uint i = 0; i < c.geometric_dimension; i++)
    coordinates[i] = 0.1*static_cast<double>(i);
  test_function f(value_size);

  // signature
  //print_scalar("signature", element.signature());

  // cell_shape
  print_scalar("cell_shape", element.cell_shape());

  // space_dimension
  print_scalar("space_dimension", element.space_dimension());

  // value_rank
  print_scalar("value_rank", element.value_rank());

  // value_dimension
  for (uint i = 0; i < element.value_rank(); i++)
    print_scalar("value_dimension", element.value_dimension(i), i);

  // evaluate_basis
  for (uint i = 0; i < element.space_dimension(); i++)
  {
    element.evaluate_basis(i, values, coordinates, c);
    print_array("evaluate_basis:", value_size, values, i);
  }

  // evaluate_basis all
  element.evaluate_basis_all(values, coordinates, c);
  print_array("evaluate_basis_all", element.space_dimension()*value_size, values);

  // evaluate_basis_derivatives
  for (uint i = 0; i < element.space_dimension(); i++)
  {
    for (uint n = 0; n <= max_derivative; n++)
    {
      uint num_derivatives = 1;
      for (uint j = 0; j < n; j++)
        num_derivatives *= c.geometric_dimension;
      element.evaluate_basis_derivatives(i, n, values, coordinates, c);
      print_array("evaluate_basis_derivatives", value_size*num_derivatives, values, i, n);
    }
  }

  // evaluate_basis_derivatives_all
  for (uint n = 0; n <= max_derivative; n++)
  {
    uint num_derivatives = 1;
      for (uint j = 0; j < n; j++)
        num_derivatives *= c.geometric_dimension;
    element.evaluate_basis_derivatives_all(n, values, coordinates, c);
    print_array("evaluate_basis_derivatives_all", element.space_dimension()*value_size*num_derivatives, values, n);
  }

  // evaluate_dof
  for (uint i = 0; i < element.space_dimension(); i++)
  {
    dof_values[i] = element.evaluate_dof(i, f, c);
    print_scalar("evaluate_dof", dof_values[i], i);
  }

  // evaluate_dofs
  element.evaluate_dofs(values, f, c);
  print_array("evaluate_dofs", element.space_dimension(), values);

  // interpolate_vertex_values
  element.interpolate_vertex_values(vertex_values, dof_values, c);
  print_array("interpolate_vertex_values", (c.topological_dimension + 1)*value_size, vertex_values);

  // num_sub_dof_elements
  print_scalar("num_sub_elements", element.num_sub_elements());

  // create_sub_element
  for (uint i = 0; i < element.num_sub_elements(); i++)
  {
    ufc::finite_element* sub_element = element.create_sub_element(i);
    test_finite_element(*sub_element);
    delete sub_element;
  }

  // Cleanup
  delete [] values;
  delete [] dof_values;
  delete [] vertex_values;
  delete [] coordinates;
}

// Function for testing ufc::element objects
void test_dofmap(ufc::dof_map& dofmap, ufc::shape cell_shape)
{
  std::cout << std::endl;
  std::cout << "Testing dof_map" << std::endl;
  std::cout << "---------------" << std::endl;

  // Prepare arguments
  test_mesh m(cell_shape);
  test_cell c(cell_shape);
  uint n = dofmap.max_local_dimension();
  uint* dofs = new uint[n];
  for (uint i = 0; i < n; i++)
    dofs[i] = 0;

  uint num_facets = c.topological_dimension + 1;
  double** coordinates = new double * [n];
  for (uint i = 0; i < n; i++)
    coordinates[i] = new double[c.geometric_dimension];

  // signature
  //print_scalar("signature", dofmap.signature());

  // needs_mesh_entities
  for (uint d = 0; d <= c.topological_dimension; d++)
    print_scalar("needs_mesh_entities", dofmap.needs_mesh_entities(d), d);

  // init_mesh
  print_scalar("init_mesh", dofmap.init_mesh(m));

  // init_cell not tested (not used by FFC)
  print_scalar("init_cell", 0);

  // init_cell_finalize (not used by FFC)
  print_scalar("init_cell_finalize", 0);

  // global_dimension
  print_scalar("global_dimension", dofmap.global_dimension());

  // local_dimension
  print_scalar("local_dimension", dofmap.local_dimension(c));

  // max_local_dimension
  print_scalar("max_local_dimension", dofmap.max_local_dimension());

  // geometric_dimension
  print_scalar("geometric_dimension", dofmap.geometric_dimension());

  // num_facet_dofs
  print_scalar("num_facet_dofs", dofmap.num_facet_dofs());

  // num_entity_dofs
  for (uint d = 0; d <= c.topological_dimension; d++)
    print_scalar("num_entity_dofs", dofmap.num_entity_dofs(d), d);

  // tabulate_dofs
  dofmap.tabulate_dofs(dofs, m, c);
  print_array("tabulate_dofs", dofmap.local_dimension(c), dofs);

  // tabulate_facet_dofs
  for (uint facet = 0; facet < num_facets; facet++)
  {
    dofmap.tabulate_facet_dofs(dofs, facet);
    print_array("tabulate_facet_dofs", dofmap.num_facet_dofs(), dofs, facet);
  }

  // tabulate_entity_dofs
  for (uint d = 0; d <= c.topological_dimension; d++)
  {
    uint num_entities[4][4] = {{0, 0, 0, 0},  // dummy entities in 0D
                               {2, 1, 0, 0},  // interval
                               {3, 3, 1, 0},  // triangle
                               {4, 6, 4, 1}}; // tetrahedron
    for (uint i = 0; i < num_entities[c.topological_dimension][d]; i++)
    {
      dofmap.tabulate_entity_dofs(dofs, d, i);
      print_array("tabulate_entity_dofs", dofmap.num_entity_dofs(d), dofs, d, i);
    }
  }

  // tabulate_coordinates
  dofmap.tabulate_coordinates(coordinates, c);
  for (uint i = 0; i < dofmap.local_dimension(c); i++)
    print_array("tabulate_coordinates", c.geometric_dimension, coordinates[i], i);

  // num_sub_dof_maps
  print_scalar("num_sub_dof_maps", dofmap.num_sub_dof_maps());

  // create_sub_dof_map
  for (uint i = 0; i < dofmap.num_sub_dof_maps(); i++)
  {
    ufc::dof_map* sub_dofmap = dofmap.create_sub_dof_map(i);
    test_dofmap(*sub_dofmap, cell_shape);
    delete sub_dofmap;
  }

  // Cleanup
  delete [] dofs;
  for (uint i = 0; i < n; i++)
    delete [] coordinates[i];
  delete [] coordinates;
}

// Function for testing ufc::cell_integral objects
void test_cell_integral(ufc::cell_integral& integral,
                        ufc::shape cell_shape,
                        int tensor_size,
                        double** w)
{
  std::cout << std::endl;
  std::cout << "Testing cell_integral" << std::endl;
  std::cout << "---------------------" << std::endl;

  // Prepare arguments
  test_cell c(cell_shape);
  double* A = new double[tensor_size];
  for(int i = 0; i < tensor_size; i++)
    A[i] = 0.0;

  // Call tabulate_tensor
  integral.tabulate_tensor(A, w, c);
  print_array("tabulate_tensor", tensor_size, A);

  // Cleanup
  delete [] A;
}

// Function for testing ufc::exterior_facet_integral objects
void test_exterior_facet_integral(ufc::exterior_facet_integral& integral,
                                  ufc::shape cell_shape,
                                  uint tensor_size,
                                  double** w)
{
  std::cout << std::endl;
  std::cout << "Testing exterior_facet_integral" << std::endl;
  std::cout << "-------------------------------" << std::endl;

  // Prepare arguments
  test_cell c(cell_shape);
  uint num_facets = c.topological_dimension + 1;
  double* A = new double[tensor_size];

  // Call tabulate_tensor for each facet
  for (uint facet = 0; facet < num_facets; facet++)
  {
    for(uint i = 0; i < tensor_size; i++)
      A[i] = 0.0;

    integral.tabulate_tensor(A, w, c, facet);
    print_array("tabulate_tensor", tensor_size, A, facet);
  }

  // Cleanup
  delete [] A;
}

// Function for testing ufc::interior_facet_integral objects
void test_interior_facet_integral(ufc::interior_facet_integral& integral,
                                  ufc::shape cell_shape,
                                  uint macro_tensor_size,
                                  double** w)
{
  std::cout << std::endl;
  std::cout << "Testing interior_facet_integral" << std::endl;
  std::cout << "-------------------------------" << std::endl;

  // Prepare arguments
  test_cell c0(cell_shape, 0);
  test_cell c1(cell_shape, 1);
  uint num_facets = c0.topological_dimension + 1;
  double* A = new double[macro_tensor_size];

  // Call tabulate_tensor for each facet-facet combination
  for (uint facet0 = 0; facet0 < num_facets; facet0++)
  {
    for (uint facet1 = 0; facet1 < num_facets; facet1++)
    {
      for(uint i = 0; i < macro_tensor_size; i++)
        A[i] = 0.0;

      integral.tabulate_tensor(A, w, c0, c1, facet0, facet1);
      print_array("tabulate_tensor", macro_tensor_size, A, facet0, facet1);
    }
  }
}

// Function for testing ufc::form objects
void test_form(ufc::form& form)
{
  std::cout << std::endl;
  std::cout << "Testing form" << std::endl;
  std::cout << "------------" << std::endl;

  // Compute size of tensors
  int tensor_size = 1;
  int macro_tensor_size = 1;
  for (uint i = 0; i < form.rank(); i++)
  {
    ufc::finite_element* element = form.create_finite_element(i);
    tensor_size *= element->space_dimension();
    macro_tensor_size *= 2*element->space_dimension();
    delete element;
  }

  // Prepare dummy coefficients
  double** w = 0;
  if (form.num_coefficients() > 0)
  {
    w = new double * [form.num_coefficients()];
    for (uint i = 0; i < form.num_coefficients(); i++)
    {
      ufc::finite_element* element = form.create_finite_element(form.rank() + i);
      w[i] = new double[element->space_dimension()];
      for (uint j = 0; j < element->space_dimension(); j++)
        w[i][j] = 0.1*static_cast<double>((i + 1)*(j + 1));
      delete element;
    }
  }

  // Get cell shape
  ufc::finite_element* element = form.create_finite_element(0);
  ufc::shape cell_shape = element->cell_shape();
  delete element;
  element = 0;

  // signature
  //print_scalar("signature", form.signature());

  // rank
  //print_scalar("rank", form.signature());

  // num_coefficients
  print_scalar("num_coefficients", form.num_coefficients());

  // num_cell_integrals
  print_scalar("num_cell_integrals", form.num_cell_integrals());

  // num_exterior_facet_integrals
  print_scalar("num_exterior_facet_integrals", form.num_exterior_facet_integrals());

  // num_interior_facet_integrals
  print_scalar("num_interior_facet_integrals", form.num_interior_facet_integrals());

  // create_finite_element
  for (uint i = 0; i < form.rank() + form.num_coefficients(); i++)
  {
    ufc::finite_element* element = form.create_finite_element(i);
    test_finite_element(*element);
    delete element;
  }

  // create_dof_map
  for (uint i = 0; i < form.rank() + form.num_coefficients(); i++)
  {
    ufc::dof_map* dofmap = form.create_dof_map(i);
    test_dofmap(*dofmap, cell_shape);
    delete dofmap;
  }

  // create_cell_integral
  for (uint i = 0; i < form.num_cell_integrals(); i++)
  {
    ufc::cell_integral* integral = form.create_cell_integral(i);
    if (integral)
      test_cell_integral(*integral, cell_shape, tensor_size, w);
    delete integral;
  }

  // create_exterior_facet_integral
  for (uint i = 0; i < form.num_exterior_facet_integrals(); i++)
  {
    ufc::exterior_facet_integral* integral = form.create_exterior_facet_integral(i);
    if (integral)
      test_exterior_facet_integral(*integral, cell_shape, tensor_size, w);
    delete integral;
  }

  // create_interior_facet_integral
  for (uint i = 0; i < form.num_interior_facet_integrals(); i++)
  {
    ufc::interior_facet_integral* integral = form.create_interior_facet_integral(i);
    if (integral)
      test_interior_facet_integral(*integral, cell_shape, macro_tensor_size, w);
    delete integral;
  }

  // Cleanup
  for (uint i = 0; i < form.num_coefficients(); i++)
    delete [] w[i];
  delete [] w;
}
