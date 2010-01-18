"""
Code snippets for testing of element/dof map code
"""

__author__ = "Marie E. Rognes (meg@simula.no)"
__license__  = "GNU GPL version 3 or any later version"


# Last modified: 15.01.2010

def element_code(element):

    ir = {}

    ir["space_dimension"] = _space_dimension()
    ir["value_dimension"] = _value_dimension()
    #ir["evaluate_basis"] = _evaluate_basis(element)
    #ir["evaluate_basis_derivatives"] = _evaluate_basis_derivatives(element)
    #ir["evaluate_dof"] = _evaluate_dof(element)

    return ir

def dofmap_code(element):

    ir = {}
    #ir["tabulate_dofs"] = _tabulate_dofs(element)
    return ir

# -------- element ------------------


template = """\
#include <iostream>
#include <ufc.h>
#include "element.h"

int main()
{

element_finite_element_0 element;

%s

return 0;
}
"""

cell_template = """
// Initialize ufc cell
ufc::cell c;
"""

function_template = """
// Initialize ufc function
ufc::function f;
"""

def _space_dimension():
    return template % """\
    int dim = element.space_dimension();
    std::cout << dim << std::endl;
    """

def _value_dimension():
    return template % """\
    int dim = element.value_dimension(0);
    std::cout << dim << std::endl;
    """

def _evaluate_basis(element):
    return template % ""

def _evaluate_basis_derivatives(element):
    return template % ""

def _evaluate_dof(element):
    " c++ code for testing evaluate_dof"

    # Declarations
    code = [cell_template, function_template]
    code += ["""
    // Initialize result:
    double result;

    // Initialize iterator
    unsigned int i;
    """]

    # Evaluate_dof code
    dofs_code = ["""\
    // Evaluate dof number %(i)d
    i = %(i)d;
    result = element.evaluate_dof(i, f, c);
    std::cout << result << std::endl;
    """ % {"i": i} for i in range(element.space_dimension())]

    code += dofs_code
    return template % "".join(code)


# -------- dof map ------------------

dof_template = """\
#include <iostream>
#include <ufc.h>
#include "element.h"

int main()
{

element_dof_map_0 dofmap;

%s

return 0;
}
"""

def _tabulate_dofs(element):
    return dof_template % """\
    //unsigned int* dofs;
    //const ufc::mesh m;

    ufc::cell cell;
    cell.coordinates = new double * [3]
    """



