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
    # ensure protocol compliance
    import numpy as np

    _formatter: interface.Formatter = Formatter(np.float64)
    # TODO: can not type compare module and class
    # _expression: interface.expression = expression
    # _file: interface.file = file
    # _integral: interface.integral = integral
