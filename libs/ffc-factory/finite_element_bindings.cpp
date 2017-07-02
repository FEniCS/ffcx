// Copyright (C) 2017 Garth N. Wells
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

#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <ufc.h>

namespace py = pybind11;

namespace ufc_wrappers
{
  void check_array_shape(py::array_t<double> x,
                         const std::vector<std::size_t> shape)
  {
    const py::buffer_info x_info = x.request();
    if ((std::size_t) x_info.ndim != shape.size())
      throw std::runtime_error("NumPy array has wrong number of dimensions");

    const auto& x_shape = x_info.shape;
    for (std::size_t i =0; i < shape.size();++i)
    {
      if ((std::size_t) x_shape[i] != shape[i])
        throw std::runtime_error("NumPy array has wrong size");
    }
  }

  // Interface function to
  // ufc::finite_element::evaluate_reference_basis
  py::array_t<double> evaluate_reference_basis(ufc::finite_element &instance,
                                               py::array_t<double> x)
  {
    // Dimensions
    const std::size_t num_dofs = instance.space_dimension();
    const std::size_t ref_size = instance.reference_value_size();
    const std::size_t tdim = instance.topological_dimension();

    // Extract point data
    const py::buffer_info x_info = x.request();
    const auto& x_shape = x_info.shape;
    const std::size_t num_points = x_shape[0];
    check_array_shape(x, {num_points, tdim});

    // Create object to hold data to be computed and returned
    py::array_t<double, py::array::c_style> f({num_points, num_dofs, ref_size});

    // Evaluate basis
    instance.evaluate_reference_basis(f.mutable_data(), num_points, x.data());

    return f;
  }

  py::array_t<double>
  evaluate_reference_basis_derivatives(const ufc::finite_element &instance,
                                       std::size_t order, py::array_t<double> x)
  {
    // Dimensions
    const std::size_t tdim = instance.topological_dimension();
    const std::size_t num_dofs = instance.space_dimension();
    const std::size_t ref_size = instance.reference_value_size();

    // Extract point data
    py::buffer_info x_info = x.request();
    const auto& x_shape = x_info.shape;
    const std::size_t num_points = x_shape[0];
    check_array_shape(x, {num_points, tdim});

    // Shape of data to be computed and returned, and create
    const std::size_t num_derivs = pow(tdim, order);
    py::array_t<double, py::array::c_style>
      f({num_points, num_dofs, num_derivs, ref_size});

    // Evaluate basis derivatives
    instance.evaluate_reference_basis_derivatives(f.mutable_data(), order,
                                                  num_points, x.data());

    return f;
  }

  py::array_t<double>
  transform_reference_basis_derivatives(ufc::finite_element &instance,
                                        std::size_t order,
                                        py::array_t<double> reference_values,
                                        py::array_t<double> X,
                                        py::array_t<double> J,
                                        py::array_t<double> detJ,
                                        py::array_t<double> K,
                                        int orientation)
  {
    // Dimensions
    const std::size_t gdim = instance.geometric_dimension();
    const std::size_t tdim = instance.topological_dimension();
    const std::size_t num_dofs = instance.space_dimension();
    const std::size_t ref_size = instance.reference_value_size();

    const std::size_t num_derivs = std::pow(tdim, order);

    // Extract reference value data (reference_values)
    py::buffer_info info_reference_values = reference_values.request();
    const auto& shape_reference_values = info_reference_values.shape;
    const std::size_t num_points = shape_reference_values[0];
    check_array_shape(reference_values, {num_points, num_dofs, num_derivs,
          ref_size});

    // Extract point data (X)
    py::buffer_info x_info = X.request();
    check_array_shape(X, {num_points, tdim});

    // Extract Jacobian data (J)
    py::buffer_info info_J = J.request();
    check_array_shape(J, {num_points, gdim, tdim});

    // Extract determinant of Jacobian data (detJ)
    py::buffer_info info_detJ = detJ.request();
    check_array_shape(detJ, {num_points});

    // Extract inverse of Jacobian data (K)
    py::buffer_info info_K = K.request();
    check_array_shape(K, {num_points, tdim, gdim});

    // Shape of data to be computed and returned, and create
    py::array_t<double, py::array::c_style>
      f({num_points, num_dofs, num_derivs, ref_size});

    // Evaluate basis derivatives
    instance.transform_reference_basis_derivatives(f.mutable_data(), order,
                                                   num_points,
                                                   reference_values.data(),
                                                   X.data(), J.data(),
                                                   detJ.data(), K.data(),
                                                   orientation);

    return f;
  }

  // Remove from UFC
  py::array_t<double> evaluate_basis(ufc::finite_element &instance,
                                     std::size_t i,
                                     py::array_t<double> x,
                                     py::array_t<double> coordinate_dofs,
                                     int cell_orientation)
  {
    // Dimensions
    const std::size_t gdim = instance.geometric_dimension();
    const std::size_t tdim = instance.topological_dimension();
    const std::size_t ref_size = instance.reference_value_size();

    // Extract point data (x)
    py::buffer_info info_x = x.request();
    check_array_shape(x, {gdim});

    // Extract coordinate data (coordinate_dofs)
    py::buffer_info info_coordinate_dofs = coordinate_dofs.request();

    // FIXME: this assumes an affine map
    check_array_shape(coordinate_dofs, {tdim + 1, gdim});

    // Shape of data to be computed and returned, and create
    py::array_t<double, py::array::c_style> values(ref_size);

    instance.evaluate_basis(i, values.mutable_data(), x.data(), coordinate_dofs.data(),
                            cell_orientation, nullptr);

    return values;
  }

  // Remove from UFC
  py::array_t<double> evaluate_basis_all(ufc::finite_element &instance,
                                         py::array_t<double> x,
                                         py::array_t<double> coordinate_dofs,
                                         int cell_orientation)
  {
    // Dimensions
    const std::size_t gdim = instance.geometric_dimension();
    const std::size_t tdim = instance.topological_dimension();
    const std::size_t num_dofs = instance.space_dimension();
    const std::size_t ref_size = instance.reference_value_size();

    // Extract point data (x)
    py::buffer_info info_x = x.request();
    check_array_shape(x, {gdim});

    // Extract reference value data (coordinate_dofs)
    py::buffer_info info_coordinate_dofs = coordinate_dofs.request();

    // FIXME: this assumes an affine map
    check_array_shape(coordinate_dofs, {tdim + 1, gdim});

    // Shape of data to be computed and returned, and create
    py::array_t<double, py::array::c_style> f({num_dofs, ref_size});

    // Call UFC function
    instance.evaluate_basis_all(f.mutable_data(), x.data(),
                                coordinate_dofs.data(),
                                cell_orientation, nullptr);

    return f;
  }

  // Remove from UFC
  py::array_t<double> evaluate_basis_derivatives(ufc::finite_element &instance,
                                                 std::size_t i,
                                                 std::size_t n,
                                                 py::array_t<double> x,
                                                 py::array_t<double> coordinate_dofs,
                                                 int cell_orientation)
  {
    // Dimensions
    const std::size_t gdim = instance.geometric_dimension();
    const std::size_t tdim = instance.topological_dimension();

    // Extract point data (x)
    py::buffer_info info_x = x.request();
    check_array_shape(x, {gdim});

    // Extract reference value data (coordinate_dofs)
    py::buffer_info info_coordinate_dofs = coordinate_dofs.request();

    // FIXME: this assumes an affine map
    check_array_shape(coordinate_dofs, {tdim + 1, gdim});

    // Shape of data to be computed and returned, and create
    py::array_t<double, py::array::c_style> f(gdim);

    // Call UFC function
    instance.evaluate_basis_derivatives(i, n, f.mutable_data(), x.data(),
                                        coordinate_dofs.data(),
                                        cell_orientation, nullptr);

    return f;
  }

  // Remove from UFC
  py::array_t<double>
  evaluate_basis_derivatives_all(ufc::finite_element &instance,
                                 std::size_t n,
                                 py::array_t<double> x,
                                 py::array_t<double> coordinate_dofs,
                                 int cell_orientation)
  {
    // Dimensions
    const std::size_t gdim = instance.geometric_dimension();
    const std::size_t tdim = instance.topological_dimension();
    const std::size_t num_dofs = instance.space_dimension();

    // Extract point data (x)
    py::buffer_info info_x = x.request();
    check_array_shape(x, {gdim});

    // Extract coordinate data (coordinate_dofs)
    py::buffer_info info_coordinate_dofs = coordinate_dofs.request();

    // FIXME: this assumes an affine map
    check_array_shape(coordinate_dofs, {tdim + 1, gdim});

    // Shape of data to be computed and returned, and create
    py::array_t<double, py::array::c_style> f({num_dofs, gdim});

    // Call UFC function
    instance.evaluate_basis_derivatives_all(n, f.mutable_data(), x.data(),
                                            coordinate_dofs.data(),
                                            cell_orientation, nullptr);

    return f;
  }

  double evaluate_dof(ufc::finite_element &instance,
                      std::size_t i,
                      const ufc::function& f,
                      py::array_t<double> coordinate_dofs,
                      int cell_orientation,
                      const ufc::cell& cell)
  {
    // Dimensions
    const std::size_t gdim = instance.geometric_dimension();
    const std::size_t tdim = instance.topological_dimension();

    // Extract coordinate data (coordinate_dofs)
    py::buffer_info info_coordinate_dofs = coordinate_dofs.request();

    // FIXME: this assumes an affine map
    check_array_shape(coordinate_dofs, {tdim + 1, gdim});

    // Call UFC function
    return instance.evaluate_dof(i, f , coordinate_dofs.data(),
                                 cell_orientation, cell, nullptr);
  }

  py::array_t<double> evaluate_dofs(ufc::finite_element &instance,
                                    const ufc::function& f,
                                    py::array_t<double> coordinate_dofs,
                                    int cell_orientation,
                                    const ufc::cell& cell)
  {
    // Dimensions
    const std::size_t gdim = instance.geometric_dimension();
    const std::size_t tdim = instance.topological_dimension();
    const std::size_t num_dofs = instance.space_dimension();

    // Extract coordinate data (coordinate_dofs)
    py::buffer_info info_coordinate_dofs = coordinate_dofs.request();

    // FIXME: this assumes an affine map
    check_array_shape(coordinate_dofs, {tdim + 1, gdim});

    // Shape of data to be computed and returned, and create
    py::array_t<double, py::array::c_style> values(num_dofs);

    // Call UFC function
    instance.evaluate_dofs(values.mutable_data(), f , coordinate_dofs.data(),
                           cell_orientation, cell, nullptr);

    return values;
  }

  py::array_t<double> interpolate_vertex_values(ufc::finite_element &instance,
                                                py::array_t<double> dof_values,
                                                py::array_t<double> coordinate_dofs,
                                                int cell_orientation)
  {
    // Dimensions
    const std::size_t gdim = instance.geometric_dimension();
    const std::size_t tdim = instance.topological_dimension();
    const std::size_t num_dofs = instance.space_dimension();
    const std::size_t num_components = instance.value_size();

    // Extract dof values (dof_values)
    py::buffer_info info_dof_values = dof_values.request();
    check_array_shape(dof_values, {num_dofs});

    // Extract coordinate data (coordinate_dofs)
    py::buffer_info info_coordinate_dofs = coordinate_dofs.request();

    // FIXME: this assumes an affine map
    check_array_shape(coordinate_dofs, {tdim + 1, gdim});

    // Shape of data to be computed and returned, and create
    std::size_t num_vertices = tdim + 1;
    py::array_t<double, py::array::c_style> vertex_values({num_vertices,
          num_components});

    // Call UFC function
    instance.interpolate_vertex_values(vertex_values.mutable_data(),
                                       dof_values.data(),
                                       coordinate_dofs.data(),
                                       cell_orientation, nullptr);

    return vertex_values;
  }

  py::array_t<double> tabulate_dof_coordinates(ufc::finite_element &instance,
                                               py::array_t<double> coordinate_dofs)
  {
    // Dimensions
    const std::size_t gdim = instance.geometric_dimension();
    const std::size_t tdim = instance.topological_dimension();
    const std::size_t num_dofs = instance.space_dimension();
    //const std::size_t num_components = instance.value_size();

    // Extract coordinate data (coordinate_dofs)
    py::buffer_info info_coordinate_dofs = coordinate_dofs.request();

    // FIXME: this assumes an affine map
    check_array_shape(coordinate_dofs, {tdim + 1, gdim});

    // Shape of data to be computed and returned, and create
    py::array_t<double, py::array::c_style> coords({num_dofs, gdim});


    // Call UFC function
    instance.tabulate_dof_coordinates(coords.mutable_data(),
                                      coordinate_dofs.data(), nullptr);

    return coords;
  }

  py::array_t<double>
  tabulate_reference_dof_coordinates(ufc::finite_element &instance)
  {
    // Dimensions
    const std::size_t gdim = instance.geometric_dimension();
    const std::size_t num_dofs = instance.space_dimension();

    // Shape of data to be computed and returned, and create
    py::array_t<double, py::array::c_style> coords({num_dofs, gdim});

    // Call UFC function
    instance.tabulate_reference_dof_coordinates(coords.mutable_data());

    return coords;
  }
}
