import os
import hashlib

# Version of FFCx header files
__author__ = "FEniCS Project"
__license__ = "This code is released into the public domain"
__version__ = "0.3.1"

# Get abspath on import, it can in some cases be a relative path w.r.t.
# curdir on startup
_include_path = os.path.dirname(os.path.abspath(__file__))


@property
def include_path():
    """Location of UFC header files"""
    return _include_path


def _compute_signature():
    # Compute signature of ufcx header file
    h = hashlib.sha1()
    with open(os.path.join(get_include_path(), "ufcx.h")) as f:
        h.update(f.read().encode("utf-8"))
    return h.hexdigest()


_signature = _compute_signature()


@property
def ufcx_sha1():
    """SHA-1 hash of the ufcx.h file

    In this implementation, the value is computed on import.
    """
    return _signature
