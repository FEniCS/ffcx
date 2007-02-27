"Code generation for cell integral"

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-02-27 -- 2007-02-27"
__copyright__ = "Copyright (C) 2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# FFC tensor representation modules
from ffc.compiler.representation.tensor import *

# FFC code generation modules
from geometrytensor import *

def generate_cell_integral(representation, format):
    """Generate dictionary of code for cell integral from the given
    form representation according to the given format"""

    code = {}

    # Generate code for tabulate_tensor
    code["tabulate_tensor"] = __generate_tabulate_tensor(representation, format)

    return code

def __generate_tabulate_tensor(representation, format):
    "Generate code for tabulate_tensor"

    # At this point, we need to check the type of representation and
    # generate code accordingly. For now, we assume that we just have
    # the tensor representation. Hint: do something differently for
    # quadrature here.

    # Generate code as a list of declarations
    code = []

    # Generate code for geometry tensor
    code += generate_geometry_tensor(representation.cell_tensor, format)

    return code
