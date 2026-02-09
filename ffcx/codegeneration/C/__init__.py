"""Generation of C code."""

from ffcx.codegeneration import interface
from ffcx.codegeneration.C import expression, file, form, integral

from .formatter import Formatter

__all__ = [
    "Formatter",
    "expression",
    "file",
    "form",
    "integral",
]

Formatter: interface.Formatter
expression: interface.expression
file: interface.file
integral: interface.integral
