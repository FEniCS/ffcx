# -*- coding: utf-8 -*-
"""
FEniCS Form Compiler (FFC)
--------------------------

FFC compiles finite element variational forms into C++ code.

The interface consists of the following functions:

  compile_form       - Compilation of forms
  compile_element    - Compilation of finite elements
  jit                - Just-In-Time compilation of forms and elements
  default_parameters - Default parameter values for FFC
  ufc_signature      - Signature of UFC interface (SHA-1 hash of ufc.h)
"""

import pkg_resources

__version__ = pkg_resources.get_distribution("fenics-ffc").version

from ffc.git_commit_hash import git_commit_hash

# Import compiler functions
from ffc.compiler import compile_form, compile_element

# Import JIT compiler
from ffc.jitcompiler import jit

# Import UFC config functions
from ffc.backends.ufc import get_include_path, get_ufc_cxx_flags, get_ufc_signature, ufc_signature

# Import default parameters
from ffc.parameters import default_parameters, default_jit_parameters

# Import plotting
from ffc.plot import *

# Duplicate list of supported elements from FIAT
from FIAT import supported_elements
supported_elements = sorted(supported_elements.keys())

# Append elements that we can plot
from ffc.plot import element_colors
supported_elements_for_plotting = list(set(supported_elements).union(set(element_colors.keys())))
supported_elements_for_plotting.sort()

# Remove elements from list that we don't support or don't trust
supported_elements.remove("Argyris")
supported_elements.remove("Hermite")
supported_elements.remove("Morley")

# Import main function, entry point to script
from ffc.main import main
