import os
import hashlib

# Version of FFC header files
__author__ = "FEniCS Project"
__license__ = "This code is released into the public domain"
__version__ = "2018.2.0.dev0"

# Get abspath on import, it can in some cases be a relative path w.r.t.
# curdir on startup
_include_path = os.path.dirname(os.path.abspath(__file__))


def get_include_path():
    """Return location of UFC header files"""
    return _include_path


def _compute_signature():
    # Compute signature of ufc header files
    h = hashlib.sha1()
    for fn in ("ufc.h", "ufc_geometry.h"):
        with open(os.path.join(get_include_path(), fn)) as f:
            h.update(f.read().encode("utf-8"))
    return h.hexdigest()


_signature = _compute_signature()


def get_signature():
    """Return SHA-1 hash of the contents of ufc.h and ufc_geometry.h.

    In this implementation, the value is computed on import.
    """
    return _signature
