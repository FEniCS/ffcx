"""Generation of numba code."""

from ffcx.codegeneration.numba import expressions, file, form, integrals

suffixes = (None, "_numba.py")

__all__ = ["expressions", "file", "form", "integrals", "suffixes"]
