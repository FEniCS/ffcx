# Start by setting the FIAT numbering scheme for entities.
# This differs between FFC and FIAT but may change in future
# versions of FIAT. It's important that we do this first
# before any other FIAT modules are loaded
try:
    from FIAT import numbering
    numbering.numbering_scheme = "UFC"
except:
    print "*** Warning: Unable to reorder entities. You need to patch or update your"
    print "*** Warning: installation of FIAT. Variable numbering_scheme is missing."
    print "*** Warning: Results may be incorrect for higher order elements on tets"

# Import finite elements and dof map
from ffc.fem.finiteelement import FiniteElement
from ffc.fem.vectorelement import VectorElement
from ffc.fem.mixedelement import MixedElement
from ffc.fem.mixedfunctions import TestFunctions, TrialFunctions, Functions
from ffc.fem.dofmap import DofMap

# Import compiler
from ffc.compiler.compiler import compile

# Import JIT compiler
from ffc.jit.jit import jit

# Import form language
from ffc.compiler.language.algebra import BasisFunction, TestFunction, TrialFunction, Function
from ffc.compiler.language.operators import *
from ffc.compiler.language.builtins import *

# Import constants
from ffc.common.constants import *
