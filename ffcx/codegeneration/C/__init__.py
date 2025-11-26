"""Generation of C code."""

from ffcx.codegeneration.C import expressions, file, form, integrals

suffixes = (".h", ".c")

__all__ = ["expressions", "file", "form", "integrals", "suffixes"]
