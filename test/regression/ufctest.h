// Copyright (C) 2010-2015 Anders Logg
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
// Functions for calling generated UFC functions with "random" (but
// fixed) data and print the output to screen. Useful for running
// regression tests.

#include <sstream>
#include <string>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <memory>
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
std::vector<double> test_coordinate_dofs(std::size_t gdim, std::size_t gdeg)
{
  // Generate some "random" coordinates
  std::vector<double> coordinate_dofs;
  if (gdim == 1)
  {
    coordinate_dofs.resize(2);
    coordinate_dofs[0]  = 0.903;
    coordinate_dofs[1]  = 0.561;
    // Huh? Only 2 vertices for interval, is this for tdim=1,gdim=2?
//    coordinate_dofs[2]  = 0.987;
//    coordinate_dofs[3]  = 0.123;
    if (gdeg > 1)
    {
      assert(gdeg == 2);
      coordinate_dofs.resize(3);
      coordinate_dofs[2]  = 0.750;
    }
  }
  else if (gdim == 2)
  {
    coordinate_dofs.resize(6);
    coordinate_dofs[0]  = 0.903;
    coordinate_dofs[1]  = 0.341;
    coordinate_dofs[2]  = 0.561;
    coordinate_dofs[3]  = 0.767;
    coordinate_dofs[4]  = 0.987;
    coordinate_dofs[5]  = 0.783;
    // Huh? Only 4 vertices for triangle, is this for quads?
//    coordinate_dofs[6]  = 0.123;
//    coordinate_dofs[7] = 0.561;
    if (gdeg > 1)
    {
      assert(gdeg == 2);
      coordinate_dofs.resize(12);
      coordinate_dofs[6]  = 0.750;
      coordinate_dofs[7]  = 0.901;
      coordinate_dofs[8]  = 0.999;
      coordinate_dofs[9]  = 0.500;
      coordinate_dofs[10] = 0.659;
      coordinate_dofs[11] = 0.555;
    }
  }
  else if (gdim == 3)
  {
    coordinate_dofs.resize(12);
    coordinate_dofs[0]  = 0.903;
    coordinate_dofs[1]  = 0.341;
    coordinate_dofs[2]  = 0.457;
    coordinate_dofs[3]  = 0.561;
    coordinate_dofs[4]  = 0.767;
    coordinate_dofs[5]  = 0.833;
    coordinate_dofs[6]  = 0.987;
    coordinate_dofs[7]  = 0.783;
    coordinate_dofs[8]  = 0.191;
    coordinate_dofs[9]  = 0.123;
    coordinate_dofs[10] = 0.561;
    coordinate_dofs[11] = 0.667;
    if (gdeg > 1)
    {
      assert(gdeg == 2);
      coordinate_dofs.resize(30);
      // FIXME: Add some quadratic tetrahedron
      assert(false);
    }
  }
  return coordinate_dofs;
}

// Function for creating "random" vertex coordinates
std::pair<std::vector<double>, std::vector<double>> test_coordinate_dof_pair(int gdim, int gdeg, int facet0, int facet1)
{
  // For affine simplices only so far...
  assert(gdeg == 1);

  // Return pair of cell coordinates there facet0 of cell 0 cooresponds to facet1 of cell 1
  int num_vertices = gdim+1;
  std::vector<double> c0 = test_coordinate_dofs(gdim, gdeg);
  std::vector<double> c1(c0.size());
  std::vector<double> m(gdim);

  for (int i=0; i<gdim; ++i)
    m[i] = 0.0;

  int k = 0;
  int j0 = 0;
  int j1 = 0;
  while (k < num_vertices-1)
  {
    // Skip the vertex that has vertex index == facet index, this is _not_ part of the facet.
    if (j0 == facet0)
      ++j0;
    if (j1 == facet1)
      ++j1;

    // Copy vertex and accumulate midpoint
    for (int i=0; i<gdim; ++i)
      m[i] += (c1[gdim*j1 + i] = c0[gdim*j0 + i]);

    // Skip to next vertex
    ++j0;
    ++j1;
    ++k;
  }
  // Scale midpoint
  for (int i=0; i<gdim; ++i)
    m[i] /= (num_vertices-1.0);

  // Place last vertex (not on facet) on the opposite side of facet
  j0 = facet0;  // Simplex specific rule
  j1 = facet1;
  for (int i=0; i<gdim; ++i)
    c1[gdim*j1 + i] = 2.0*m[i] - c0[gdim*j0 + i];

  return {c0, c1};
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
    case ufc::shape::interval:
      topological_dimension = 1;
      if (gdim == 0)
        geometric_dimension = 1;
      else
        geometric_dimension = gdim;
      break;
    case ufc::shape::triangle:
      topological_dimension = 2;
      if (gdim == 0)
        geometric_dimension = 2;
      else
        geometric_dimension = gdim;
      break;
    case ufc::shape::tetrahedron:
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
    entity_indices.resize(topological_dimension+1);
    for (std::size_t i = 0; i <= topological_dimension; i++)
    {
      entity_indices[i].resize(12); // hack: sufficient for any cell entities list, 12 edges of hex
      for (std::size_t j = 0; j < entity_indices[i].size(); j++)
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
  // NOTE: Assuming geometry degree 1
  const std::vector<double> coordinate_dofs
    = test_coordinate_dofs(element.geometric_dimension(), 1);
  std::size_t value_size = 1;
  for (std::size_t i = 0; i < element.value_rank(); i++)
    value_size *= element.value_dimension(i);
  std::size_t derivative_size = 1;
  for (std::size_t i = 0; i < max_derivative; i++)
    derivative_size *= c.geometric_dimension;
  std::vector<double>
    values(element.space_dimension()*value_size*derivative_size, 0.0);

  std::vector<double> dof_values(element.space_dimension(), 0.0);
  std::vector<double> vertex_values((c.topological_dimension + 1)*value_size,
                                    0.0);

  std::vector<double> coordinates(c.geometric_dimension);
  for (std::size_t i = 0; i < c.geometric_dimension; i++)
    coordinates[i] = 0.1*static_cast<double>(i);
  test_function f(value_size);

  // signature
  //printer.print_scalar("signature", element.signature());

  // cell_shape
  printer.print_scalar("cell_shape", static_cast<int>(element.cell_shape()));

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
    element.evaluate_basis(i, values.data(), coordinates.data(),
                           coordinate_dofs.data(), 1);
    printer.print_array("evaluate_basis:", value_size, values.data(), i);
  }

// evaluate_basis all
  element.evaluate_basis_all(values.data(), coordinates.data(),
                             coordinate_dofs.data(), 1);
  printer.print_vector("evaluate_basis_all", values);

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
                                         values.data(),
                                         coordinates.data(),
                                         coordinate_dofs.data(),
                                         1);
      printer.print_array("evaluate_basis_derivatives",
                          value_size*num_derivatives, values.data(), i, n);
    }
  }

  // evaluate_basis_derivatives_all
  for (std::size_t n = 0; n <= max_derivative; n++)
  {
    std::size_t num_derivatives = 1;
    for (std::size_t j = 0; j < n; j++)
      num_derivatives *= c.geometric_dimension;
    element.evaluate_basis_derivatives_all(n,
                                           values.data(),
                                           coordinates.data(),
                                           coordinate_dofs.data(),
                                           1);
    printer.print_array("evaluate_basis_derivatives_all",
                        element.space_dimension()*value_size*num_derivatives,
                        values.data(), n);
  }

  // evaluate_dof
  for (std::size_t i = 0; i < element.space_dimension(); i++)
  {
    dof_values[i] = element.evaluate_dof(i, f,
                                         coordinate_dofs.data(),
                                         1, c);
    printer.print_scalar("evaluate_dof", dof_values[i], i);
  }

// evaluate_dofs
  element.evaluate_dofs(values.data(), f, coordinate_dofs.data(), 1, c);
  printer.print_array("evaluate_dofs", element.space_dimension(),
                       values.data());

  // interpolate_vertex_values
  element.interpolate_vertex_values(vertex_values.data(),
                                    dof_values.data(),
                                    coordinate_dofs.data(),
                                    1);
  printer.print_vector("interpolate_vertex_values", vertex_values);

  // tabulate_coordinates
  std::vector<double>
    dof_coordinates(element.geometric_dimension()*element.space_dimension());
  element.tabulate_dof_coordinates(dof_coordinates.data(),
                                   coordinate_dofs.data());
  printer.print_vector("tabulate_dof_coordinates", dof_coordinates);

  // num_sub_dof_elements
  printer.print_scalar("num_sub_elements", element.num_sub_elements());

  // create_sub_element
  for (std::size_t i = 0; i < element.num_sub_elements(); i++)
  {
    std::unique_ptr<ufc::finite_element>
      sub_element(element.create_sub_element(i));
    test_finite_element(*sub_element, i, printer);
  }

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

  //test_cell c(cell_shape, dofmap.geometric_dimension());
  test_cell c(cell_shape, dofmap.topological_dimension());
  std::size_t n = dofmap.num_element_dofs();
  std::vector<std::size_t> dofs(n, 0);

  std::size_t num_facets = c.topological_dimension + 1;
  std::vector<double> coordinates(n*c.geometric_dimension);

  // needs_mesh_entities
  for (std::size_t d = 0; d <= c.topological_dimension; d++)
  {
    printer.print_scalar("needs_mesh_entities", dofmap.needs_mesh_entities(d),
                         d);
  }

  // global_dimension
  printer.print_scalar("global_dimension",
                       dofmap.global_dimension(num_entities));

  // num_element_dofs
  printer.print_scalar("num_element_dofs", dofmap.num_element_dofs());

  // num_facet_dofs
  printer.print_scalar("num_facet_dofs", dofmap.num_facet_dofs());

  // num_entity_dofs
  for (std::size_t d = 0; d <= c.topological_dimension; d++)
    printer.print_scalar("num_entity_dofs", dofmap.num_entity_dofs(d), d);

  // tabulate_dofs
  dofmap.tabulate_dofs(dofs.data(), num_entities, c.entity_indices);
  printer.print_array("tabulate_dofs", dofmap.num_element_dofs(), dofs.data());

  // tabulate_facet_dofs
  for (std::size_t facet = 0; facet < num_facets; facet++)
  {
    dofmap.tabulate_facet_dofs(dofs.data(), facet);
    printer.print_array("tabulate_facet_dofs", dofmap.num_facet_dofs(),
                        dofs.data(), facet);
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
      dofmap.tabulate_entity_dofs(dofs.data(), d, i);
      printer.print_array("tabulate_entity_dofs", dofmap.num_entity_dofs(d),
                          dofs.data(), d, i);
    }
  }

  // num_sub_dofmaps
  printer.print_scalar("num_sub_dofmaps", dofmap.num_sub_dofmaps());

  // create_sub_dofmap
  for (std::size_t i = 0; i < dofmap.num_sub_dofmaps(); i++)
  {
    std::unique_ptr<ufc::dofmap> sub_dofmap(dofmap.create_sub_dofmap(i));
    test_dofmap(*sub_dofmap, cell_shape, i, printer);
  }

  printer.end();
}

// Function for testing ufc::cell_integral objects
void test_cell_integral(ufc::cell_integral& integral,
                        ufc::shape cell_shape,
                        std::size_t gdim,
                        std::size_t gdeg,
                        std::size_t tensor_size,
                        double** w,
                        bool bench,
                        int id,
                        Printer & printer)
{
  printer.begin("cell_integral", id);

  // Prepare arguments
  test_cell c(cell_shape, gdim);
  const std::vector<double> coordinate_dofs = test_coordinate_dofs(gdim, gdeg);
  std::vector<double> A(tensor_size, 0.0);

  // Call tabulate_tensor
  integral.tabulate_tensor(A.data(), w, coordinate_dofs.data(), c.orientation);
  printer.print_vector("tabulate_tensor", A);

  // Benchmark tabulate tensor
  if (bench)
  {
    printer.begin("timing");
    for (std::size_t num_reps = initial_num_reps;; num_reps *= 2)
    {
      double t0 = time();
      for (std::size_t i = 0; i < num_reps; i++)
      {
        integral.tabulate_tensor(A.data(), w, coordinate_dofs.data(),
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

  printer.end();
}

// Function for testing ufc::exterior_facet_integral objects
void test_exterior_facet_integral(ufc::exterior_facet_integral& integral,
                                  ufc::shape cell_shape,
                                  std::size_t gdim,
                                  std::size_t gdeg,
                                  std::size_t tensor_size,
                                  double** w,
                                  bool bench,
                                  int id,
                                  Printer & printer)
{
  printer.begin("exterior_facet_integral", id);

  // Prepare arguments
  test_cell c(cell_shape, gdim);
  const std::vector<double> coordinate_dofs = test_coordinate_dofs(gdim, gdeg);
  std::size_t num_facets = c.topological_dimension + 1;
  std::vector<double> A(tensor_size);

  // Call tabulate_tensor for each facet
  for (std::size_t facet = 0; facet < num_facets; facet++)
  {
    for(std::size_t i = 0; i < tensor_size; i++)
      A[i] = 0.0;

    integral.tabulate_tensor(A.data(), w, coordinate_dofs.data(), facet,
                             c.orientation);
    printer.print_array("tabulate_tensor", tensor_size, A.data(), facet);
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
        integral.tabulate_tensor(A.data(), w, coordinate_dofs.data(), 0,
                                 c.orientation);
      }
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

  printer.end();
}

// Function for testing ufc::interior_facet_integral objects
void test_interior_facet_integral(ufc::interior_facet_integral& integral,
                                  ufc::shape cell_shape,
                                  std::size_t gdim,
                                  std::size_t gdeg,
                                  std::size_t macro_tensor_size,
                                  double** w,
                                  bool bench,
                                  int id,
                                  Printer & printer)
{
  printer.begin("interior_facet_integral", id);

  // Prepare arguments
  test_cell c0(cell_shape, gdim);
  //test_cell c1(cell_shape, gdim);
  std::size_t num_facets = c0.topological_dimension + 1;
  std::vector<double> A(macro_tensor_size);

  // Call tabulate_tensor for each facet-facet combination
  for (std::size_t facet0 = 0; facet0 < num_facets; facet0++)
  {
    for (std::size_t facet1 = 0; facet1 < num_facets; facet1++)
    {
      const std::pair<std::vector<double>, std::vector<double>> coordinate_dofs
        = test_coordinate_dof_pair(gdim, gdeg, facet0, facet1);

      for(std::size_t i = 0; i < macro_tensor_size; i++)
        A[i] = 0.0;

      integral.tabulate_tensor(A.data(),
                               w,
                               coordinate_dofs.first.data(),
                               coordinate_dofs.second.data(),
                               facet0, facet1,
                               1, 1);
      printer.print_array("tabulate_tensor", macro_tensor_size, A.data(),
                          facet0, facet1);
    }
  }

  // Benchmark tabulate tensor
  if (bench)
  {
    const std::pair<std::vector<double>, std::vector<double>> coordinate_dofs
      = test_coordinate_dof_pair(gdim, gdeg, 0, 0);

    printer.begin("timing");
    for (std::size_t num_reps = initial_num_reps;; num_reps *= 2)
    {
      double t0 = time();
      for (std::size_t i = 0; i < num_reps; i++)
      {
        integral.tabulate_tensor(A.data(), w,
                                 coordinate_dofs.first.data(),
                                 coordinate_dofs.second.data(),
                                 0, 0,
                                 1, 1);
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

  printer.end();
}

// Function for testing ufc::vertex_integral objects
void test_vertex_integral(ufc::vertex_integral& integral,
                         ufc::shape cell_shape,
                         std::size_t gdim,
                         std::size_t gdeg,
                         std::size_t tensor_size,
                         double** w,
                         bool bench,
                         int id,
                         Printer & printer)
{
  printer.begin("vertex_integral", id);

  // Prepare arguments
  test_cell c(cell_shape, gdim);
  const std::vector<double> coordinate_dofs = test_coordinate_dofs(gdim, gdeg);
  std::size_t num_vertices = c.topological_dimension + 1;
  std::vector<double> A(tensor_size);

  // Call tabulate_tensor for each vertex
  for (std::size_t vertex = 0; vertex < num_vertices; vertex++)
  {
    for(std::size_t i = 0; i < tensor_size; i++)
      A[i] = 0.0;

    integral.tabulate_tensor(A.data(), w, coordinate_dofs.data(), vertex,
                             c.orientation);
    printer.print_array("tabulate_tensor", tensor_size, A.data(), vertex);
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
        integral.tabulate_tensor(A.data(), w, coordinate_dofs.data(), 0,
                                 c.orientation);
      }
      double dt = time() - t0;
      if (dt > minimum_timing)
      {
        dt /= static_cast<double>(num_reps);
        printer.print_scalar("vertex_integral_timing_iterations", num_reps);
        printer.print_scalar("vertex_integral_time", dt);
        break;
      }
    }
    printer.end();
  }

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
    std::unique_ptr<ufc::finite_element> element(form.create_finite_element(i));
    tensor_size *= element->space_dimension();

    // *2 for interior facet integrals
    macro_tensor_size *= 2*element->space_dimension();
  }

  // Prepare dummy coefficients
  double** w = 0;
  if (form.num_coefficients() > 0)
  {
    w = new double * [form.num_coefficients()];
    for (std::size_t i = 0; i < form.num_coefficients(); i++)
    {
      std::unique_ptr<ufc::finite_element>
        element(form.create_finite_element(form.rank() + i));

      // *2 for interior facet integrals
      const std::size_t macro_dim = 2*element->space_dimension();

      w[i] = new double[macro_dim];
      for (std::size_t j = 0; j < macro_dim; j++)
        w[i][j] = 0.1*static_cast<double>((i + 1)*(j + 1));
    }
  }

  // Get cell shape
  std::unique_ptr<ufc::finite_element> element(form.create_coordinate_finite_element());
  ufc::shape cell_shape = element->cell_shape();
  std::size_t gdim = element->geometric_dimension();
  std::size_t gdeg = element->degree();
  assert(element->value_rank() == 1);
  assert(element->value_dimension(0) == gdim);
  assert(element->family() == std::string("Lagrange"));
  element.reset();

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

  // has_vertex_integrals
  printer.print_scalar("has_vertex_integrals", form.has_vertex_integrals());

  // max_cell_subdomain_id
  printer.print_scalar("max_cell_subdomain_id", form.max_cell_subdomain_id());

  // max_exterior_facet_subdomain_id
  printer.print_scalar("max_exterior_facet_subdomain_id",
                       form.max_exterior_facet_subdomain_id());

  // max_interior_facet_subdomain_id
  printer.print_scalar("max_interior_facet_subdomain_id",
                       form.max_interior_facet_subdomain_id());

  // max_vertex_subdomain_id
  printer.print_scalar("max_vertex_subdomain_id",
                       form.max_vertex_subdomain_id());

  // create_finite_element
  for (std::size_t i = 0; i < form.rank() + form.num_coefficients(); i++)
  {
    std::unique_ptr<ufc::finite_element> element(form.create_finite_element(i));
    test_finite_element(*element, i, printer);
  }

  // create_dofmap
  for (std::size_t i = 0; i < form.rank() + form.num_coefficients(); i++)
  {
    std::unique_ptr<ufc::dofmap> dofmap(form.create_dofmap(i));
    test_dofmap(*dofmap, cell_shape, i, printer);
  }

  // create_cell_integral
  {
    std::unique_ptr<ufc::cell_integral>
      integral(form.create_default_cell_integral());
    printer.print_scalar("default_cell_integral", (bool)integral);
    if (integral)
    {
      test_cell_integral(*integral, cell_shape, gdim, gdeg,
                         tensor_size, w, bench, -1, printer);
    }
  }
  for (std::size_t i = 0; i < form.max_cell_subdomain_id(); i++)
  {
    std::unique_ptr<ufc::cell_integral> integral(form.create_cell_integral(i));
    if (integral)
    {
      test_cell_integral(*integral, cell_shape, gdim, gdeg,
                         tensor_size, w, bench, i, printer);
    }
  }

  // create_exterior_facet_integral
  {
    std::unique_ptr<ufc::exterior_facet_integral>
      integral(form.create_default_exterior_facet_integral());
    printer.print_scalar("default_exterior_facet_integral", (bool)integral);
    if (integral)
    {
      test_exterior_facet_integral(*integral, cell_shape, gdim, gdeg,
                                   tensor_size, w, bench, -1, printer);
    }
  }

  for (std::size_t i = 0; i < form.max_exterior_facet_subdomain_id(); i++)
  {
    std::unique_ptr<ufc::exterior_facet_integral>
      integral(form.create_exterior_facet_integral(i));
    if (integral)
    {
      test_exterior_facet_integral(*integral, cell_shape, gdim, gdeg,
                                   tensor_size, w, bench, i, printer);
    }
  }

  // create_interior_facet_integral
  {
    std::unique_ptr<ufc::interior_facet_integral>
      integral(form.create_default_interior_facet_integral());
    printer.print_scalar("default_interior_facet_integral", (bool)integral);
    if (integral)
    {
      test_interior_facet_integral(*integral, cell_shape, gdim, gdeg,
                                   macro_tensor_size, w, bench, -1, printer);
    }
  }
  for (std::size_t i = 0; i < form.max_interior_facet_subdomain_id(); i++)
  {
    std::unique_ptr<ufc::interior_facet_integral>
      integral(form.create_interior_facet_integral(i));
    if (integral)
    {
      test_interior_facet_integral(*integral, cell_shape, gdim, gdeg,
                                   macro_tensor_size, w, bench, i, printer);
    }
  }

  // create_vertex_integral
  {
    std::unique_ptr<ufc::vertex_integral>
      integral(form.create_default_vertex_integral());
    printer.print_scalar("default_vertex_integral", (bool)integral);
    if (integral)
    {
      test_vertex_integral(*integral, cell_shape, gdim, gdeg,
                           tensor_size, w, bench, -1, printer);
    }
  }
  for (std::size_t i = 0; i < form.max_vertex_subdomain_id(); i++)
  {
    std::unique_ptr<ufc::vertex_integral>
      integral(form.create_vertex_integral(i));
    if (integral)
    {
      test_vertex_integral(*integral, cell_shape, gdim, gdeg,
                           tensor_size, w, bench, i, printer);
    }
  }

  // Cleanup
  for (std::size_t i = 0; i < form.num_coefficients(); i++)
    delete [] w[i];
  delete [] w;

  printer.end();
}
