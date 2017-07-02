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

#include <memory>
#include <pybind11/pybind11.h>
#include <ufc.h>

namespace py = pybind11;

namespace ufc_wrappers
{
  void finite_element(py::module& m);
  void dofmap(py::module& m);
}

//PYBIND11_MODULE(ffc_factory, m)
PYBIND11_PLUGIN(ffc_factory)
{
  pybind11::module m("ffc_factory", "Factory for wrapping FFC JIT-compiled objects");

  // Create finite_element submodule
  py::module finite_element = m.def_submodule("finite_element",
                                              "UFC finite_element");
  ufc_wrappers::finite_element(finite_element);

  // Create dofmap submodule
  py::module dofmap = m.def_submodule("dofmap", "UFC dofmap");
  ufc_wrappers::dofmap(dofmap);

  // Factory functions
  m.def("make_ufc_finite_element",
        [](std::uintptr_t e)
        {
          ufc::finite_element * p = reinterpret_cast<ufc::finite_element *>(e);
          return std::shared_ptr<const ufc::finite_element>(p);
        });
  m.def("make_ufc_dofmap",
        [](std::uintptr_t e)
        {
          ufc::dofmap * p = reinterpret_cast<ufc::dofmap *>(e);
          return std::shared_ptr<const ufc::dofmap>(p);
        });

  return m.ptr();
}
