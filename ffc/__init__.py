# Start by setting the FIAT numbering scheme for entities and the
# reference element. This differs between FFC and FIAT but may change
# in future versions of FIAT. It's important that we do this first
# before any other FIAT modules are loaded

from FIAT import numbering
numbering.numbering_scheme = "UFC"

from FIAT import reference
reference.reference_element = "default"

# Import compiler
from ffc.compiler.uflcompiler import compile as compile_ufl

# Import JIT compiler
from ffc.jit.jit import jit
