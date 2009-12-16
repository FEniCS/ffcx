"""
FEniCS Form Compiler (FFC)
--------------------------

FFC compiles finite element variational forms into C++ code.

The interface consists of the following three functions:

  compile_form    - Compilation of forms
  compile_element - Compilation of finite elements
  jit             - Just-In-Time compilation of forms and elements
"""

# Import compiler functions
from ffc.compiler import compile_form, compile_element

# Import JIT compiler
from ffc.jit import jit
