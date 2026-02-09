"""Generation of numba code."""

from typing import TYPE_CHECKING

from ffcx.codegeneration import interface
from ffcx.codegeneration.numba import expression, file, form, integral

from .formatter import Formatter

__all__ = [
    "Formatter",
    "expression",
    "file",
    "form",
    "integral",
]

if TYPE_CHECKING:
    Formatter: interface.Formatter
    expression: interface.expression
    file: interface.file
    integral: interface.integral
