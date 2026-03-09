"""Generation of C code."""

from typing import TYPE_CHECKING

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

if TYPE_CHECKING:
    # ensure protocol compliance
    import numpy as np

    _formatter: interface.Formatter = Formatter(np.float64)
    _expression: interface.expression_generator = expression.generator
    _file: interface.file_generator = file.generator
    _form: interface.form_generator = form.generator
    _integral: interface.integral_generator = integral.generator
