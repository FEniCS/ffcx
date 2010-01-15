"""
Code snippets for for generated element and dof map code
"""

__author__ = "Marie E. Rognes (meg@simula.no)"
__license__  = "GNU GPL version 3 or any later version"


# Last modified: 15.01.2010

element_common = """\
#include <iostream>
#include "element.h"

int main()
{

element_finite_element_0 element;

%s

return 0;
}
"""

space_dimension = """\
int dim = element.space_dimension();
std::cout << dim << std::endl;
"""

value_dimension = """\
int dim = element.value_dimension(0);
std::cout << dim << std::endl;
"""

space_dimension = element_common % space_dimension
value_dimension = element_common % value_dimension
