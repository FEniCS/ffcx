# -*- coding: utf-8 -*-
# Based on original implementation by Martin Alnes and Anders Logg
#
# Modified by Anders Logg 2015

__all__ = ["dolfin_tag", "stl_includes", "dolfin_includes", "snippets"]

dolfin_tag = "// DOLFIN wrappers"

stl_includes = """\
// Standard library includes
#include <string>
"""

dolfin_includes = """\
// DOLFIN includes
#include <dolfin/common/NoDeleter.h>
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/MultiMesh.h>
#include <dolfin/fem/FiniteElement.h>
#include <dolfin/fem/DofMap.h>
#include <dolfin/fem/Form.h>
#include <dolfin/fem/MultiMeshForm.h>
#include <dolfin/function/FunctionSpace.h>
#include <dolfin/function/MultiMeshFunctionSpace.h>
#include <dolfin/function/GenericFunction.h>
#include <dolfin/function/CoefficientAssigner.h>
#include <dolfin/function/MultiMeshCoefficientAssigner.h>
#include <dolfin/adaptivity/ErrorControl.h>
#include <dolfin/adaptivity/GoalFunctional.h>
#include <dolfin/la/GenericVector.h>"""


snippets = {"shared_ptr_space":
            ("std::shared_ptr<const dolfin::FunctionSpace> %s",
             "    _function_spaces[%d] = %s;"),
            "referenced_space":
            ("const dolfin::FunctionSpace& %s",
             "    _function_spaces[%d] = reference_to_no_delete_pointer(%s);"),
            "multimesh_shared_ptr_space":
            ("std::shared_ptr<const dolfin::MultiMeshFunctionSpace> %s",
             None),
            "multimesh_referenced_space":
            ("const dolfin::MultiMeshFunctionSpace& %s",
             None),
            "shared_ptr_mesh":
            ("std::shared_ptr<const dolfin::Mesh> mesh",
             "    _mesh = mesh;"),
            "shared_ptr_multimesh":
                ("std::shared_ptr<const dolfin::MultiMesh> mesh",
                 "    _multimesh = mesh;"),
            "referenced_mesh":
            ("const dolfin::Mesh& mesh",
             "    _mesh = reference_to_no_delete_pointer(mesh);"),
            "shared_ptr_coefficient":
            ("std::shared_ptr<const dolfin::GenericFunction> %s",
             "    this->%s = %s;"),
            "shared_ptr_ref_coefficient":
            ("std::shared_ptr<const dolfin::GenericFunction> %s",
             "    this->%s = *%s;"),
            "referenced_coefficient":
            ("const dolfin::GenericFunction& %s",
             "    this->%s = %s;"),
            "functionspace":
            ("TestSpace", "TrialSpace"),
            "multimeshfunctionspace":
            ("MultiMeshTestSpace", "MultiMeshTrialSpace")
            }
