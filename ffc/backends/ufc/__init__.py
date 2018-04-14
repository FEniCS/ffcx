# -*- coding: utf-8 -*-
"""Utility functions for UFC"""

__author__ = "FEniCS Project"
__version__ = "2018.1.0.dev0"
__license__ = "This code is released into the public domain"

import os
from hashlib import sha1

# Get abspath on import, it can in some cases be a relative path w.r.t.
# curdir on startup
_include_path = os.path.dirname(os.path.abspath(__file__))


def get_include_path():
    "Return location of UFC header files"
    return _include_path


def _compute_ufc_signature():
    # Compute signature of ufc header files
    h = sha1()
    for fn in ("ufc.h", "ufc_geometry.h"):
        with open(os.path.join(get_include_path(), fn)) as f:
            h.update(f.read().encode("utf-8"))
    return h.hexdigest()


_ufc_signature = _compute_ufc_signature()


def get_ufc_signature():
    """Return SHA-1 hash of the contents of ufc.h and ufc_geometry.h.

    In this implementation, the value is computed on import.
    """
    return _ufc_signature


def get_ufc_cxx_flags():
    """Return C++ flags for compiling UFC C++11 code.

    Return type is a list of strings.

    Used internally in some tests.
    """
    return ["-std=c++14"]


# ufc_signature() already introduced to FFC standard in 1.7.0dev,
# called by the dolfin cmake build system to compare against
# future imported ffc versions for compatibility.
ufc_signature = get_ufc_signature
