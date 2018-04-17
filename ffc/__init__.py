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

import logging

import pkg_resources

# Import JIT compiler
from ffc.jitcompiler import jit  # noqa: F401

# Import main function, entry point to script
from ffc.main import main  # noqa: F401

# Import default parameters
from ffc.parameters import (default_jit_parameters, default_parameters)  # noqa: F401
from FIAT import supported_elements

__version__ = pkg_resources.get_distribution("fenics-ffc").version

logging.basicConfig()
logger = logging.getLogger("ffc")


class FFCError(Exception):
    pass
    # def __init__(self,*args,**kwargs):
    #       raise FFCError("LLLLLLLLLL", self)
    #       Exception.__init__(self, *args, **kwargs)


# from ffc.git_commit_hash import git_commit_hash

# Import compiler functions
# from ffc.compiler import compile_form, compile_element

# Duplicate list of supported elements from FIAT
supported_elements = sorted(supported_elements.keys())

# Remove elements from list that we don't support or don't trust
supported_elements.remove("Argyris")
supported_elements.remove("Hermite")
supported_elements.remove("Morley")
