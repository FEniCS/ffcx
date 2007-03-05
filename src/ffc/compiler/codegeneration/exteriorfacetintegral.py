"Code generation for exterior facet integral"

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2007-02-28 -- 2007-03-05"
__copyright__ = "Copyright (C) 2007 Anders Logg"
__license__  = "GNU GPL Version 2"

# FFC tensor representation modules
from ffc.compiler.representation.tensor import *

# FFC code generation modules
import tensorrepresentation

def generate_exterior_facet_integral(representation, sub_domain, format):
    """Generate dictionary of code for exterior facet integral from
    the given form representation according to the given format"""

    code = {}

    # Generate code for tabulate_tensor
    code["tabulate_tensor"] = __generate_tabulate_tensor(representation, sub_domain, format)

    return code

def __generate_tabulate_tensor(representation, sub_domain, format):
    "Generate code for tabulate_tensor"

    # At this point, we need to check the type of representation and
    # generate code accordingly. For now, we assume that we just have
    # the tensor representation. Hint: do something differently for
    # quadrature here.

    return tensorrepresentation.generate_tabulate_tensor(representation, sub_domain, format)
