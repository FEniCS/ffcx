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
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#include <ufc.h>

namespace py = pybind11;

namespace ufc_wrappers
{

  void dofmap(py::module& m)
  {
    // Wrap ufc::dofmap classe
    py::class_<ufc::dofmap, std::shared_ptr<ufc::dofmap>>
      dofmap(m, "dofmap");

    dofmap.def("topological_dimension",
               &ufc::dofmap::topological_dimension);
  }

}
