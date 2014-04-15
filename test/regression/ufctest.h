// Copyright (C) 2010-2013 Anders Logg
//
// This file is part of FFC.
//
// FFC is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// FFC is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with FFC. If not, see <http://www.gnu.org/licenses/>.
//
// Modified by Martin Alnaes, 2013
//
// First added:  2010-01-24
// Last changed: 2013-02-18
//
// Functions for calling generated UFC functions with "random" (but
// fixed) data and print the output to screen. Useful for running
// regression tests.

#include <sstream>
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <ufc.h>

#include "printer.h"

// How many derivatives to test
const std::size_t max_derivative = 2;

// Parameters for adaptive timing
const std::size_t initial_num_reps = 10;
const double minimum_timing = 1.0;

// Function for timing
double time()
{
  clock_t __toc_time = std::clock();
  return ((double) (__toc_time)) / CLOCKS_PER_SEC;
}

// Function for creating "random" vertex coordinates
std::vector<double> test_vertex_coordinates(int gdim)
{
  // Generate some "random" coordinates
  std::vector<double> vertex_coordinates;
  if (gdim == 1)
  {
    vertex_coordinates.resize(4);
    vertex_coordinates[0]  = 0.903;
    vertex_coordinates[1]  = 0.561;
    vertex_coordinates[2]  = 0.987;
    vertex_coordinates[3]  = 0.123;
  }
  else if (gdim == 2)
  {
    vertex_coordinates.resize(8);
    vertex_coordinates[0]  = 0.903;
    vertex_coordinates[1]  = 0.341;
    vertex_coordinates[2]  = 0.561;
    vertex_coordinates[3]  = 0.767;
    vertex_coordinates[4]  = 0.987;
    vertex_coordinates[5]  = 0.783;
    vertex_coordinates[6]  = 0.123;
    vertex_coordinates[7] = 0.561;
  }
  else if (gdim == 3)
  {
    vertex_coordinates.resize(12);
    vertex_coordinates[0]  = 0.903;
    vertex_coordinates[1]  = 0.341;
    vertex_coordinates[2]  = 0.457;
    vertex_coordinates[3]  = 0.561;
    vertex_coordinates[4]  = 0.767;
    vertex_coordinates[5]  = 0.833;
    vertex_coordinates[6]  = 0.987;
    vertex_coordinates[7]  = 0.783;
    vertex_coordinates[8]  = 0.191;
    vertex_coordinates[9]  = 0.123;
    vertex_coordinates[10] = 0.561;
    vertex_coordinates[11] = 0.667;
  }
  return vertex_coordinates;
}

// Class for creating "random" ufc::cell objects
class test_cell : public ufc::cell
{
public:

  test_cell(ufc::shape cell_shape, int gdim=0, int offset=0)
  {
    // Store cell shape
    this->cell_shape = cell_shape;

    // Store dimensions
    switch (cell_shape)
    {
    case ufc::interval:
      topological_dimension = 1;
      if (gdim == 0)
        geometric_dimension = 1;
      else
        geometric_dimension = gdim;
      break;
    case ufc::triangle:
      topological_dimension = 2;
      if (gdim == 0)
        geometric_dimension = 2;
      else
        geometric_dimension = gdim;
      break;
    case ufc::tetrahedron:
      topological_dimension = 3;
      if (gdim == 0)
        geometric_dimension = 3;
      else
        geometric_dimension = gdim;
      break;
    default:
      throw std::runtime_error("Unhandled cell shape.");
    }

    // Set orientation (random, but must be set)
    this->orientation = 1;

    // Generate some "random" entity indices
    entity_indices.resize(4);
    for (std::size_t i = 0; i < 4; i++)
    {
      entity_indices[i].resize(6);
      for (std::size_t j = 0; j < 6; j++)
        entity_indices[i][j] = i*j + offset;
    }
  }

  ~test_cell() {}

};

// Class for creating a "random" ufc::function object
class test_function : public ufc::function
{
public:

  test_function(std::size_t value_size) : value_size(value_size) {}

  void evaluate(double* values, const double* coordinates,
                const ufc::cell& c) const
  {
    for (std::size_t i = 0; i < value_size; i++)
    {
      values[i] = 1.0;
      for (std::size_t j = 0; j < c.geometric_dimension; j++)
        values[i] *= static_cast<double>(i + 1)*coordinates[j];
    }
  }

private:

  std::size_t value_size;

};

std::string format_name(std::string name, int i=-1, int j=-1)
{
  std::stringstream s;
  s << name;
  if (i >= 0)
    s << "_" << i;
  if (j >= 0)
    s << "_" << j;
  return s.str();
}

// Function for testing ufc::element objects
void test_finite_element(ufc::finite_element& element, int id, Printer& printer)
{
  printer.begin("finite_element", id);

  // Prepare arguments
  test_cell c(element.cell_shape(), element.geometric_dimension());
  const std::vector<double> vertex_coordinates
    = test_vertex_coordinates(element.geometric_dimension());
  std::size_t value_size = 1;
  for (std::size_t i = 0; i < element.value_rank(); i++)
    value_size *= element.value_dimension(i);
  std::size_t derivative_size = 1;
  for (std::size_t i = 0; i < max_derivative; i++)
    derivative_size *= c.geometric_dimension;
  double* values
    = new double[element.space_dimension()*value_size*derivative_size];
  for (std::size_t i = 0;
       i < element.space_dimension()*value_size*derivative_size; i++)
  {
    values[i] = 0.0;
  }
  double* dof_values = new double[element.space_dimension()];
  for (std::size_t i = 0; i < element.space_dimension(); i++)
    dof_values[i] = 0.0;

  double* vertex_values = new double[(c.topological_dimension + 1)*value_size];
  for (std::size_t i = 0; i < (c.topological_dimension + 1)*value_size; i++)
    vertex_values[i] = 0.0;

  double* coordinates = new double[c.geometric_dimension];
  for (std::size_t i = 0; i < c.geometric_dimension; i++)
    coordinates[i] = 0.1*static_cast<double>(i);
  test_function f(value_size);

  // signature
  //printer.print_scalar("signature", element.signature());

  // cell_shape
  printer.print_scalar("cell_shape", element.cell_shape());

  // space_dimension
  printer.print_scalar("space_dimension", element.space_dimension());

  // value_rank
  printer.print_scalar("value_rank", element.value_rank());

  // value_dimension
  for (std::size_t i = 0; i < element.value_rank(); i++)
    printer.print_scalar("value_dimension", element.value_dimension(i), i);

  // evaluate_basis
  for (std::size_t i = 0; i < element.space_dimension(); i++)
  {
    element.evaluate_basis(i, values, coordinates,
                           vertex_coordinates.data(), 1);
    printer.print_array("evaluate_basis:", value_size, values, i);
  }

  // evaluate_basis all
  element.evaluate_basis_all(values, coordinates,
                             vertex_coordinates.data(), 1);
  printer.print_array("evaluate_basis_all",
                      element.space_dimension()*value_size, values);

  // evaluate_basis_derivatives
  for (std::size_t i = 0; i < element.space_dimension(); i++)
  {
    for (std::size_t n = 0; n <= max_derivative; n++)
    {
      std::size_t num_derivatives = 1;
      for (std::size_t j = 0; j < n; j++)
        num_derivatives *= c.geometric_dimension;
      element.evaluate_basis_derivatives(i,
                                         n,
                                         values,
                                         coordinates,
                                         vertex_coordinates.data(),
                                         1);
      printer.print_array("evaluate_basis_derivatives",
                          value_size*num_derivatives, values, i, n);
    }
  }

  // evaluate_basis_derivatives_all
  for (std::size_t n = 0; n <= max_derivative; n++)
  {
    std::size_t num_derivatives = 1;
    for (std::size_t j = 0; j < n; j++)
      num_derivatives *= c.geometric_dimension;
    element.evaluate_basis_derivatives_all(n,
                                           values,
                                           coordinates,
                                           vertex_coordinates.data(),
                                           1);
    printer.print_array("evaluate_basis_derivatives_all",
                        element.space_dimension()*value_size*num_derivatives,
                        values, n);
  }

  // evaluate_dof
  for (std::size_t i = 0; i < element.space_dimension(); i++)
  {
    dof_values[i] = element.evaluate_dof(i, f,
                                         vertex_coordinates.data(),
                                         1, c);
    printer.print_scalar("evaluate_dof", dof_values[i], i);
  }

  // evaluate_dofs
  element.evaluate_dofs(values, f, vertex_coordinates.data(), 1, c);
  printer.print_array("evaluate_dofs", element.space_dimension(), values);

  // interpolate_vertex_values
  element.interpolate_vertex_values(vertex_values,
                                    dof_values,
                                    vertex_coordinates.data(),
                                    1,
                                    c);
  printer.print_array("interpolate_vertex_values",
                      (c.topological_dimension + 1)*value_size,
                      vertex_values);

  // num_sub_dof_elements
  printer.print_scalar("num_sub_elements", element.num_sub_elements());

  // create_sub_element
  for (std::size_t i = 0; i < element.num_sub_elements(); i++)
  {
    ufc::finite_element* sub_element = element.create_sub_element(i);
    test_finite_element(*sub_element, i, printer);
    delete sub_element;
  }

  // Cleanup
  delete [] values;
  delete [] dof_values;
  delete [] vertex_values;
  delete [] coordinates;

  printer.end();
}

// Function for testing ufc::element objects
void test_dofmap(ufc::dofmap& dofmap, ufc::shape cell_shape, int id,
                 Printer& printer)
{
  printer.begin("dofmap", id);

  // Prepare arguments
  std::vector<std::size_t> num_entities(4);
  num_entities[0] = 10001;
  num_entities[1] = 10002;
  num_entities[2] = 10003;
  num_entities[3] = 10004;

  test_cell c(cell_shape, dofmap.geometric_dimension());
  const std::vector<double> vertex_coordinates
    = test_vertex_coordinates(dofmap.geometric_dimension());
  std::size_t n = dofmap.local_dimension();
  std::size_t* dofs = new std::size_t[n];
  for (std::size_t i = 0; i < n; i++)
    dofs[i] = 0;

  std::size_t num_facets = c.topological_dimension + 1;
  double** coordinates = new double * [n];
  for (std::size_t i = 0; i < n; i++)
    coordinates[i] = new double[c.geometric_dimension];

  // signature
  //printer.print_scalar("signature", dofmap.signature());

  // needs_mesh_entities
  for (std::size_t d = 0; d <= c.topological_dimension; d++)
  {
    printer.print_scalar("needs_mesh_entities", dofmap.needs_mesh_entities(d),
                         d);
  }

  // global_dimension
  printer.print_scalar("global_dimension",
                       dofmap.global_dimension(num_entities));

  // local_dimension
  printer.print_scalar("local_dimension", dofmap.local_dimension());

  // geometric_dimension
  printer.print_scalar("geometric_dimension", dofmap.geometric_dimension());

  // num_facet_dofs
  printer.print_scalar("num_facet_dofs", dofmap.num_facet_dofs());

  // num_entity_dofs
  for (std::size_t d = 0; d <= c.topological_dimension; d++)
    printer.print_scalar("num_entity_dofs", dofmap.num_entity_dofs(d), d);

  // tabulate_dofs
  dofmap.tabulate_dofs(dofs, num_entities, c);
  printer.print_array("tabulate_dofs", dofmap.local_dimension(), dofs);

  // tabulate_facet_dofs
  for (std::size_t facet = 0; facet < num_facets; facet++)
  {
    dofmap.tabulate_facet_dofs(dofs, facet);
    printer.print_array("tabulate_facet_dofs", dofmap.num_facet_dofs(), dofs,
                        facet);
  }

  // tabulate_entity_dofs
  for (std::size_t d = 0; d <= c.topological_dimension; d++)
  {
    std::size_t num_entities[4][4] = {{0, 0, 0, 0},  // dummy entities in 0D
                               {2, 1, 0, 0},  // interval
                               {3, 3, 1, 0},  // triangle
                               {4, 6, 4, 1}}; // tetrahedron
    for (std::size_t i = 0; i < num_entities[c.topological_dimension][d]; i++)
    {
      dofmap.tabulate_entity_dofs(dofs, d, i);
      printer.print_array("tabulate_entity_dofs", dofmap.num_entity_dofs(d),
                          dofs, d, i);
    }
  }

  // tabulate_coordinates
  dofmap.tabulate_coordinates(coordinates, vertex_coordinates.data());
  for (std::size_t i = 0; i < dofmap.local_dimension(); i++)
  {
    printer.print_array("tabulate_coordinates", c.geometric_dimension,
                        coordinates[i], i);
  }

  // num_sub_dofmaps
  printer.print_scalar("num_sub_dofmaps", dofmap.num_sub_dofmaps());

  // create_sub_dofmap
  for (std::size_t i = 0; i < dofmap.num_sub_dofmaps(); i++)
  {
    ufc::dofmap* sub_dofmap = dofmap.create_sub_dofmap(i);
    test_dofmap(*sub_dofmap, cell_shape, i, printer);
    delete sub_dofmap;
  }

  // Cleanup
  delete [] dofs;
  for (std::size_t i = 0; i < n; i++)
    delete [] coordinates[i];
  delete [] coordinates;

  printer.end();
}

// Function for testing ufc::cell_integral objects
void test_cell_integral(ufc::cell_integral& integral,
                        ufc::shape cell_shape,
                        std::size_t gdim,
                        std::size_t tensor_size,
                        double** w,
                        bool bench,
                        int id,
                        Printer & printer)
{
  printer.begin("cell_integral", id);

  // Prepare arguments
  test_cell c(cell_shape, gdim);
  const std::vector<double> vertex_coordinates = test_vertex_coordinates(gdim);
  double* A = new double[tensor_size];
  for(std::size_t i = 0; i < tensor_size; i++)
    A[i] = 0.0;

  // Call tabulate_tensor
  integral.tabulate_tensor(A, w, vertex_coordinates.data(), c.orientation);
  printer.print_array("tabulate_tensor", tensor_size, A);

  // Benchmark tabulate tensor
  if (bench)
  {
    printer.begin("timing");
    for (std::size_t num_reps = initial_num_reps;; num_reps *= 2)
    {
      double t0 = time();
      for (std::size_t i = 0; i < num_reps; i++)
      {
        integral.tabulate_tensor(A, w, vertex_coordinates.data(),
                                 c.orientation);
      }
      double dt = time() - t0;
      if (dt > minimum_timing)
      {
        dt /= static_cast<double>(num_reps);
        printer.print_scalar("cell_integral_timing_iterations", num_reps);
        printer.print_scalar("cell_integral_time", dt);
        break;
      }
    }
    printer.end();
  }

  // Cleanup
  delete [] A;

  printer.end();
}

// Function for testing ufc::exterior_facet_integral objects
void test_exterior_facet_integral(ufc::exterior_facet_integral& integral,
                                  ufc::shape cell_shape,
                                  std::size_t gdim,
                                  std::size_t tensor_size,
                                  double** w,
                                  bool bench,
                                  int id,
                                  Printer & printer)
{
  printer.begin("exterior_facet_integral", id);

  // Prepare arguments
  test_cell c(cell_shape, gdim);
  const std::vector<double> vertex_coordinates = test_vertex_coordinates(gdim);
  std::size_t num_facets = c.topological_dimension + 1;
  double* A = new double[tensor_size];

  // Call tabulate_tensor for each facet
  for (std::size_t facet = 0; facet < num_facets; facet++)
  {
    for(std::size_t i = 0; i < tensor_size; i++)
      A[i] = 0.0;

    integral.tabulate_tensor(A, w, vertex_coordinates.data(), facet,
                             c.orientation);
    printer.print_array("tabulate_tensor", tensor_size, A, facet);
  }

  // Benchmark tabulate tensor
  if (bench)
  {
    printer.begin("timing");
    for (std::size_t num_reps = initial_num_reps;; num_reps *= 2)
    {
      double t0 = time();
      for (std::size_t i = 0; i < num_reps; i++)
        integral.tabulate_tensor(A, w, vertex_coordinates.data(), 0,
                                 c.orientation);
      double dt = time() - t0;
      if (dt > minimum_timing)
      {
        dt /= static_cast<double>(num_reps);
        printer.print_scalar("exterior_facet_integral_timing_iterations",
                             num_reps);
        printer.print_scalar("exterior_facet_integral_time", dt);
        break;
      }
    }
    printer.end();
  }

  // Cleanup
  delete [] A;

  printer.end();
}

// Function for testing ufc::interior_facet_integral objects
void test_interior_facet_integral(ufc::interior_facet_integral& integral,
                                  ufc::shape cell_shape,
                                  std::size_t gdim,
                                  std::size_t macro_tensor_size,
                                  double** w,
                                  bool bench,
                                  int id,
                                  Printer & printer)
{
  printer.begin("interior_facet_integral", id);

  // Prepare arguments
  test_cell c0(cell_shape, gdim, 0);
  test_cell c1(cell_shape, gdim, 1);
  const std::vector<double> vertex_coordinates0
    = test_vertex_coordinates(gdim);
  const std::vector<double> vertex_coordinates1
    = test_vertex_coordinates(gdim);
  std::size_t num_facets = c0.topological_dimension + 1;
  double* A = new double[macro_tensor_size];

  // Call tabulate_tensor for each facet-facet combination
  for (std::size_t facet0 = 0; facet0 < num_facets; facet0++)
  {
    for (std::size_t facet1 = 0; facet1 < num_facets; facet1++)
    {
      for(std::size_t i = 0; i < macro_tensor_size; i++)
        A[i] = 0.0;

      integral.tabulate_tensor(A,
                               w,
                               vertex_coordinates0.data(),
                               vertex_coordinates1.data(),
                               facet0, facet1,
                               c0.orientation,
                               c1.orientation);
      printer.print_array("tabulate_tensor", macro_tensor_size, A,
                          facet0, facet1);
    }
  }

  // Benchmark tabulate tensor
  if (bench)
  {
    printer.begin("timing");
    for (std::size_t num_reps = initial_num_reps;; num_reps *= 2)
    {
      double t0 = time();
      for (std::size_t i = 0; i < num_reps; i++)
      {
        integral.tabulate_tensor(A, w,
                                 vertex_coordinates0.data(),
                                 vertex_coordinates1.data(), 0, 0,
                                 c0.orientation,
                                 c1.orientation);
      }

      double dt = time() - t0;
      if (dt > minimum_timing)
      {
        dt /= static_cast<double>(num_reps);
        printer.print_scalar("interior_facet_integral_timing_iterations",
                             num_reps);
        printer.print_scalar("interior_facet_integral_time", dt);
        break;
      }
    }
    printer.end();
  }

  // Cleanup
  delete [] A;

  printer.end();
}

// Function for testing ufc::point_integral objects
void test_point_integral(ufc::point_integral& integral,
                         ufc::shape cell_shape,
                         std::size_t gdim,
                         std::size_t tensor_size,
                         double** w,
                         bool bench,
                         int id,
                         Printer & printer)
{
  printer.begin("point_integral", id);

  // Prepare arguments
  test_cell c(cell_shape, gdim);
  const std::vector<double> vertex_coordinates = test_vertex_coordinates(gdim);
  std::size_t num_vertices = c.topological_dimension + 1;
  double* A = new double[tensor_size];

  // Call tabulate_tensor for each vertex
  for (std::size_t vertex = 0; vertex < num_vertices; vertex++)
  {
    for(std::size_t i = 0; i < tensor_size; i++)
      A[i] = 0.0;

    integral.tabulate_tensor(A, w, vertex_coordinates.data(), vertex,
                             c.orientation);
    printer.print_array("tabulate_tensor", tensor_size, A, vertex);
  }

  // Benchmark tabulate tensor
  if (bench)
  {
    printer.begin("timing");
    for (std::size_t num_reps = initial_num_reps;; num_reps *= 2)
    {
      double t0 = time();
      for (std::size_t i = 0; i < num_reps; i++)
        integral.tabulate_tensor(A, w, vertex_coordinates.data(), 0,
                                 c.orientation);
      double dt = time() - t0;
      if (dt > minimum_timing)
      {
        dt /= static_cast<double>(num_reps);
        printer.print_scalar("point_integral_timing_iterations", num_reps);
        printer.print_scalar("point_integral_time", dt);
        break;
      }
    }
    printer.end();
  }

  // Cleanup
  delete [] A;

  printer.end();
}

// Function for testing ufc::form objects
void test_form(ufc::form& form, bool bench, int id, Printer & printer)
{
  printer.begin("form", id);

  // Compute size of tensors
  int tensor_size = 1;
  int macro_tensor_size = 1;
  for (std::size_t i = 0; i < form.rank(); i++)
  {
    ufc::finite_element* element = form.create_finite_element(i);
    tensor_size *= element->space_dimension();
    macro_tensor_size *= 2*element->space_dimension(); // *2 for interior facet integrals
    delete element;
  }

  // Prepare dummy coefficients
  double** w = 0;
  if (form.num_coefficients() > 0)
  {
    w = new double * [form.num_coefficients()];
    for (std::size_t i = 0; i < form.num_coefficients(); i++)
    {
      ufc::finite_element* element
        = form.create_finite_element(form.rank() + i);
      const std::size_t macro_dim = 2*element->space_dimension(); // *2 for interior facet integrals
      w[i] = new double[macro_dim];
      for (std::size_t j = 0; j < macro_dim; j++)
        w[i][j] = 0.1*static_cast<double>((i + 1)*(j + 1));
      delete element;
    }
  }

  // Get cell shape
  ufc::finite_element* element = form.create_finite_element(0);
  ufc::shape cell_shape = element->cell_shape();
  std::size_t gdim = element->geometric_dimension();
  delete element;
  element = 0;

  // signature
  //printer.print_scalar("signature", form.signature());

  // rank
  printer.print_scalar("rank", form.rank());

  // num_coefficients
  printer.print_scalar("num_coefficients", form.num_coefficients());

  // has_cell_integrals
  printer.print_scalar("has_cell_integrals", form.has_cell_integrals());

  // has_exterior_facet_integrals
  printer.print_scalar("has_exterior_facet_integrals",
                       form.has_exterior_facet_integrals());

  // has_interior_facet_integrals
  printer.print_scalar("has_interior_facet_integrals",
                       form.has_interior_facet_integrals());

  // has_point_integrals
  printer.print_scalar("has_point_integrals", form.has_point_integrals());

  // num_cell_domains
  printer.print_scalar("num_cell_domains", form.num_cell_domains());

  // num_exterior_facet_domains
  printer.print_scalar("num_exterior_facet_domains",
                       form.num_exterior_facet_domains());

  // num_interior_facet_domains
  printer.print_scalar("num_interior_facet_domains",
                       form.num_interior_facet_domains());

  // num_point_domains
  printer.print_scalar("num_point_domains", form.num_point_domains());

  // create_finite_element
  for (std::size_t i = 0; i < form.rank() + form.num_coefficients(); i++)
  {
    ufc::finite_element* element = form.create_finite_element(i);
    test_finite_element(*element, i, printer);
    delete element;
  }

  // create_dofmap
  for (std::size_t i = 0; i < form.rank() + form.num_coefficients(); i++)
  {
    ufc::dofmap* dofmap = form.create_dofmap(i);
    test_dofmap(*dofmap, cell_shape, i, printer);
    delete dofmap;
  }

  // create_cell_integral
  {
    ufc::cell_integral* integral = form.create_default_cell_integral();
    printer.print_scalar("default_cell_integral", (bool)integral);
    if (integral)
      test_cell_integral(*integral, cell_shape, gdim,
                         tensor_size, w, bench, -1, printer);
    delete integral;
  }
  for (std::size_t i = 0; i < form.num_cell_domains(); i++)
  {
    ufc::cell_integral* integral = form.create_cell_integral(i);
    if (integral)
      test_cell_integral(*integral, cell_shape, gdim,
                         tensor_size, w, bench, i, printer);
    delete integral;
  }

  // create_exterior_facet_integral
  {
    ufc::exterior_facet_integral* integral
      = form.create_default_exterior_facet_integral();
    printer.print_scalar("default_exterior_facet_integral", (bool)integral);
    if (integral)
    {
      test_exterior_facet_integral(*integral, cell_shape, gdim,
                                   tensor_size, w, bench, -1, printer);
    }
    delete integral;
  }

  for (std::size_t i = 0; i < form.num_exterior_facet_domains(); i++)
  {
    ufc::exterior_facet_integral* integral
      = form.create_exterior_facet_integral(i);
    if (integral)
    {
      test_exterior_facet_integral(*integral, cell_shape, gdim,
                                   tensor_size, w, bench, i, printer);
    }
    delete integral;
  }

  // create_interior_facet_integral
  {
    ufc::interior_facet_integral* integral = form.create_default_interior_facet_integral();
    printer.print_scalar("default_interior_facet_integral", (bool)integral);
    if (integral)
      test_interior_facet_integral(*integral, cell_shape, gdim,
                                   macro_tensor_size, w, bench, -1, printer);
    delete integral;
  }
  for (std::size_t i = 0; i < form.num_interior_facet_domains(); i++)
  {
    ufc::interior_facet_integral* integral
      = form.create_interior_facet_integral(i);
    if (integral)
    {
      test_interior_facet_integral(*integral, cell_shape, gdim,
                                   macro_tensor_size, w, bench, i, printer);
    }
    delete integral;
  }

  // create_point_integral
  {
    ufc::point_integral* integral = form.create_default_point_integral();
    printer.print_scalar("default_point_integral", (bool)integral);
    if (integral)
    {
      test_point_integral(*integral, cell_shape, gdim, tensor_size, w, bench,
                          -1, printer);
    }
    delete integral;
  }
  for (std::size_t i = 0; i < form.num_point_domains(); i++)
  {
    ufc::point_integral* integral = form.create_point_integral(i);
    if (integral)
    {
      test_point_integral(*integral, cell_shape, gdim, tensor_size, w, bench,
                          i, printer);
    }
    delete integral;
  }

  // Cleanup
  for (std::size_t i = 0; i < form.num_coefficients(); i++)
    delete [] w[i];
  delete [] w;

  printer.end();
}
