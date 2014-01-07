"""
FEniCS Form Compiler (FFC)
--------------------------

FFC compiles finite element variational forms into C++ code.

The interface consists of the following functions:

  compile_form       - Compilation of forms
  compile_element    - Compilation of finite elements
  jit                - Just-In-Time compilation of forms and elements
  default_parameters - Default parameter values for FFC
"""

__version__ = "1.3.0"

# Import compiler functions
from ffc.compiler import compile_form, compile_element

# Import JIT compiler
from ffc.jitcompiler import jit

# Import default parameters
from parameters import default_parameters

# Import plotting
from plot import *

# Import useful extra functionality
from extras import *

# List of supported elements
try:

    # Import list of supported elements from FIAT
    from FIAT import supported_elements
    supported_elements = supported_elements.keys()
    supported_elements.sort()

    # Append elements that we can plot
    from plot import element_colors
    supported_elements_for_plotting = list(set(supported_elements).union(set(element_colors.keys())))
    supported_elements_for_plotting.sort()

    # Remove elements from list that we don't support or don't trust
    supported_elements.remove("Argyris")
    supported_elements.remove("Hermite")
    supported_elements.remove("Morley")

except:

    supported_elements = []
    supported_elements_for_plotting = []
