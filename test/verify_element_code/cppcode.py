"""
Code snippets for testing of element/dof map code
"""

__author__ = "Marie E. Rognes (meg@simula.no)"
__license__  = "GNU GPL version 3 or any later version"

import numpy

# Last modified: 15.01.2010

def test_code(element):

    test = {}

    # Element stuff
    test["space_dimension"] = _space_dimension()
    test["value_dimension"] = _value_dimension()
    #test["evaluate_basis"] = _evaluate_basis(element)
    #test["evaluate_basis_derivatives"] = _evaluate_basis_derivatives(element)
    test["evaluate_dof"] = _evaluate_dof(element)
    test["interpolate_vertex_values"] = _interpolate_vertex_values(element)
    #test["num_sub_elements"] = _num_sub_elements(element)

    # Dofmap stuff
    test["num_facet_dofs"] = _num_facet_dofs(element)
    #test["num_entity_dofs"] = _num_entity_dofs(element)
    test["tabulate_dofs"] = _tabulate_dofs(element)
    test["tabulate_facet_dofs"] = _tabulate_facet_dofs(element)
    test["tabulate_coordinates"] = _tabulate_coordinates(element)
    return test


# --------- general utilities ---------------

def extract_value_dimension(element):
    try:
        N = element.value_shape()[0]
    except:
        N = 1
    return N

def extract_cell_dimension(element):
    try:
        return element.cell().geometric_dimension()
    except:
        return element.ref_el.get_spatial_dimension()


def _mesh(d):
    if d == 1:
        return "dolfin::UnitInterval m(5);\nm.init();\ndolfin::UFCMesh mesh(m);"
    elif d == 2:
        return "dolfin::UnitSquare m(5, 7);\nm.init();\ndolfin::UFCMesh mesh(m);"
    elif d == 3:
        return "dolfin::UnitCube m(5, 6, 3);\nm.init();\ndolfin::UFCMesh mesh(m);"

def _cell(i=None):
    i = i or 3
    return "dolfin::UFCCell cell(dolfin::Cell(m, %s));" % i

def _function(n):

    values = ["values[%d] = sin(%d*x[0]);" % (i, i+1) for i in range(n)]
    if n == 1:
        value_declaration = ""
    else:
        value_declaration = "public:\n  Source() : Expression(%d) {}\n\n" % n

    code = """
class Source : public dolfin::Expression
{
  %s

  void eval(dolfin::Array<double>& values, const dolfin::Array<const double>& x) const
  {
    %s
   }
};

Source f;
""" % (value_declaration, "\n".join(values))
    return code

dolfin_includes = """
#include <dolfin/common/Array.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/UnitInterval.h>
#include <dolfin/mesh/UnitSquare.h>
#include <dolfin/mesh/UnitCube.h>
#include <dolfin/function/Expression.h>
#include <dolfin/fem/UFC.h>
"""

# -------- element ------------------

template = """\
#include <iostream>
%s
#include "element.h"

int main()
{

element_0_finite_element_0 element;

%s

return 0;
}
"""

def _space_dimension():
    "Test code for space_dimension."
    code = """\
int dim = element.space_dimension();
std::cout << dim << std::endl;"""
    return template % ("", code)

def _value_dimension():
    "Test code for value_dimension."

    code = """\
int dim = element.value_dimension(0);
std::cout << dim << std::endl;"""
    return template % ("", code)

def _evaluate_basis(element):
    "Test code for evaluate_basis."

    return template % ("", code)

def _evaluate_basis_derivatives(element):
    "Test code for evaluate_basis_derivatives."

    return template % ("", code)

def _evaluate_dof(element):
    "Test code for evaluate_dof."

    # Initialize cell
    d = extract_cell_dimension(element)
    code = [_mesh(d), _cell()]

    # Initialize function f
    n = extract_value_dimension(element)
    code += [_function(n)]

    # Evaluate_dof code
    code += ["double result;"]
    for i in range(element.space_dimension()):
        code += ["result = element.evaluate_dof(%d, f, cell);" % i]
        code += ["printf(\"%0.5g, \", result);"]

    return template % (dolfin_includes, "\n".join(code))


def _interpolate_vertex_values(element):
    "Test code for interpolate_vertex_values."

    N = extract_value_dimension(element)
    d = extract_cell_dimension(element)
    n = element.space_dimension()

    # Initialize default mesh and cell
    code = [_mesh(d), _cell()]

    # Initialize dof values
    dof_values = ["%d" % i for i in range(n)]
    code += ["double dof_values[%d] = {%s};" % (n, ",".join(dof_values))]

    # Declare vertex_values:
    code += ["double vertex_values[%d];" % (N*(d+1))]

    # Call interpolate
    func = "element.interpolate_vertex_values"
    code += [func + "(vertex_values, dof_values, cell);"]

    # Print resulting vertex values
    for i in range(N*(d+1)):
        code += ["printf(\"%0.5g, \"" + ", vertex_values[%d]);" % i]

    return template % (dolfin_includes, "\n".join(code))

def _num_sub_elements(element):
    "Test code for num_sub_elements."

    code = ["int num_elements = element.num_sub_elements();"]
    code += ["std::cout << num_elements << std::endl;"]
    return template % ("", "\n".join(code))

# -------- dof map ------------------

dof_template = """\
#include <iostream>
%s
#include "element.h"

int main()
{

element_0_dof_map_0 dofmap;

%s

return 0;
}
"""


def _num_facet_dofs(element):
    "Test code for num_facet_dofs."

    code = "unsigned int result = dofmap.num_facet_dofs();"
    return dof_template % ("", code)

def _num_entity_dofs(element):
    "Test code for num_entity_dofs."

    d = extract_cell_dimension(element)

    code = ["unsigned int result;"]
    for i in range(d):
        code += ["result = dofmap.num_entity_dofs(%d);" % i]
        code += ["std::cout << result << std::endl;"]
    return dof_template % ("", "\n".join(code))


def _tabulate_dofs(element):
    "Test code for tabulate_dofs."

    d = extract_cell_dimension(element)
    code = []

    # Initialize default mesh and cell
    code += [_mesh(d), _cell()]

    # Declare dof array
    n = element.space_dimension()
    code += ["unsigned int dofs[%d];" % n]

    # Tabulate dofs
    code += ["dofmap.tabulate_dofs(dofs, mesh, cell);"]

    # Print dofs
    code += ["std::cout << dofs[%d] << std::endl;" % i for i in range(n)]

    return dof_template % (dolfin_includes, "\n".join(code))

def _tabulate_coordinates(element):
    "Test code for tabulate_coordinates."

    n = element.space_dimension()
    d = extract_cell_dimension(element)
    code = []

    # Initialize default mesh and cell
    code += [_mesh(d), _cell()]

    # Initialize coords
    code += ["double** coords = new double*[%d];" % n]
    for i in range(n):
        code += ["coords[%d] = new double[%d];" % (i, d)]

    # Tabulate coordinates and print
    code += ["dofmap.tabulate_coordinates(coords, cell);"]
    for i in range(n):
        for j in range(d):
            code += ["printf(\"%0.5g, \"" +", coords[%d][%d]);" % (i, j)]

    return dof_template % (dolfin_includes, "\n".join(code))

def _tabulate_facet_dofs(element):
    "Test code for tabulate_facet_dofs."

    code = []

    # Map (spatial dimension - 1) to number of facets
    facets = [2, 3, 6]

    # Declare dofs
    code += ["int n = dofmap.num_facet_dofs();"]
    code += ["unsigned int* dofs = new unsigned int[n];"]

    D = extract_cell_dimension(element) - 1
    for d in range(facets[D]):
        code += ["""
        dofmap.tabulate_facet_dofs(dofs, %d);
        for (unsigned int i=0; i < n; i++) {
          std::cout << dofs[i] << std::endl;
        }
        """ % d]

    return dof_template % ("", "\n".join(code))
