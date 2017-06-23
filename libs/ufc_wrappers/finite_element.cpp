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
  // Interface function to
  // ufc::finite_element::evaluate_reference_basis
  py::array_t<double> evaluate_reference_basis(ufc::finite_element &instance,
                                               py::array_t<double> x);

  py::array_t<double>
  evaluate_reference_basis_derivatives(const ufc::finite_element &instance,
                                       std::size_t order, py::array_t<double> x);

  py::array_t<double>
  transform_reference_basis_derivatives(ufc::finite_element &instance,
                                        std::size_t order,
                                        py::array_t<double> reference_values,
                                        py::array_t<double> X,
                                        py::array_t<double> J,
                                        py::array_t<double> detJ,
                                        py::array_t<double> K,
                                        int orientation);

  // Remove from UFC
  py::array_t<double>
  evaluate_basis(ufc::finite_element &instance, std::size_t i,
                 py::array_t<double> x,
                 py::array_t<double> coordinate_dofs,
                 int cell_orientation);

  // Remove from UFC
  py::array_t<double> evaluate_basis_all(ufc::finite_element &instance,
                                         py::array_t<double> x,
                                         py::array_t<double> coordinate_dofs,
                                         int cell_orientation);
  // Remove from UFC
  py::array_t<double> evaluate_basis_derivatives(ufc::finite_element &instance,
                                                 std::size_t i,
                                                 std::size_t n,
                                                 py::array_t<double> x,
                                                 py::array_t<double> coordinate_dofs,
                                                 int cell_orientation);

  // Remove from UFC
  py::array_t<double>
  evaluate_basis_derivatives_all(ufc::finite_element &instance,
                                 std::size_t n,
                                 py::array_t<double> x,
                                 py::array_t<double> coordinate_dofs,
                                 int cell_orientation);

  double evaluate_dof(ufc::finite_element &instance,
                      std::size_t i,
                      const ufc::function& f,
                      py::array_t<double> coordinate_dofs,
                      int cell_orientation,
                      const ufc::cell& cell);

  py::array_t<double> evaluate_dofs(ufc::finite_element &instance,
                                    const ufc::function& f,
                                    py::array_t<double> coordinate_dofs,
                                    int cell_orientation,
                                    const ufc::cell& cell);

  py::array_t<double> interpolate_vertex_values(ufc::finite_element &instance,
                                                py::array_t<double> dof_values,
                                                py::array_t<double> coordinate_dofs,
                                                int cell_orientation);

  py::array_t<double> tabulate_dof_coordinates(ufc::finite_element &instance,
                                               py::array_t<double> coordinate_dofs);

  py::array_t<double>
  tabulate_reference_dof_coordinates(ufc::finite_element &instance);

  // Add the Python submodule 'finite_element'
  void finite_element(py::module& m)
  {
    // Wrap ufc classes
    py::class_<ufc::finite_element, std::shared_ptr<ufc::finite_element>>
      element(m, "finite_element");

    // Wrap ufc::shape enum
    py::enum_<ufc::shape> shape(m, "shape");
    shape.value("interval", ufc::shape::interval);
    shape.value("triangle", ufc::shape::triangle);
    shape.value("tetrahedon", ufc::shape::tetrahedron);
    shape.value("quadrilateral", ufc::shape::quadrilateral);
    shape.value("hexahedron", ufc::shape::hexahedron);
    shape.value("vertex", ufc::shape::vertex);

    // Expose ufc::finite_element functions
    element.def("signature", &ufc::finite_element::signature);
    element.def("cell_shape", &ufc::finite_element::cell_shape);
    element.def("topological_dimension", &ufc::finite_element::topological_dimension);
    element.def("geometric_dimension", &ufc::finite_element::geometric_dimension);
    element.def("space_dimension", &ufc::finite_element::space_dimension);
    element.def("value_rank", &ufc::finite_element::value_rank);
    element.def("value_dimension", &ufc::finite_element::value_dimension);
    element.def("value_size", &ufc::finite_element::value_size);
    element.def("reference_value_rank", &ufc::finite_element::reference_value_rank);
    element.def("reference_value_dimension", &ufc::finite_element::reference_value_dimension);
    element.def("reference_value_size", &ufc::finite_element::reference_value_size);
    element.def("degree", &ufc::finite_element::degree);
    element.def("family", &ufc::finite_element::family);

    element.def("evaluate_reference_basis", &ufc_wrappers::evaluate_reference_basis,
                "evaluate reference basis");
    element.def("evaluate_reference_basis_derivatives",
                &ufc_wrappers::evaluate_reference_basis_derivatives,
              "evaluate reference basis");
    element.def("transform_reference_basis_derivatives",
                &ufc_wrappers::transform_reference_basis_derivatives,
                "transform reference basis derivatives");

    element.def("evaluate_basis", &ufc_wrappers::evaluate_basis, "evaluate basis");
    element.def("evaluate_basis_all", &ufc_wrappers::evaluate_basis_all,
                "evaluate basis all");

    element.def("evaluate_basis_derivatives", &ufc_wrappers::evaluate_basis_derivatives,
                "evaluate basis derivatives");
    element.def("evaluate_basis_derivatives_all", &ufc_wrappers::evaluate_basis_derivatives_all,
                "evaluate basis derivatives all");

    element.def("evaluate_dof", &ufc_wrappers::evaluate_dof, "evaluate dof");
    element.def("evaluate_dofs", &ufc_wrappers::evaluate_dofs, "evaluate dofs");

    element.def("interpolate_vertex_values", &ufc_wrappers::interpolate_vertex_values,
                "interpolate vertex values");
    element.def("tabulate_dof_coordinates", &ufc_wrappers::tabulate_dof_coordinates,
                "tabulate dof coordinates");

    element.def("tabulate_reference_dof_coordinates",
                &ufc_wrappers::tabulate_reference_dof_coordinates,
                "tabulate reference dof coordinates");

    element.def("num_sub_elements", &ufc::finite_element::num_sub_elements);

    element.def("create_sub_element", &ufc::finite_element::create_sub_element,
                py::return_value_policy::take_ownership);
    element.def("create", &ufc::finite_element::create,
                py::return_value_policy::take_ownership);
  }

}
