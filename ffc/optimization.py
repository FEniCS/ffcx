"""
Compiler stage 5: optimization
------------------------------

This module implements the optimization of an intermediate code
representation.
"""

__author__ = "Anders Logg (logg@simula.no)"
__date__ = "2009-12-22"
__copyright__ = "Copyright (C) 2009 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Last changed: 2009-12-22

# FFC modules
from ffc.log import info, begin, end

def optimize_ir(ir, form_data):
    "Optimize intermediate form representation."

    begin("Compiler stage 3: Optimizing intermediate representation")
    info("Optimization is currently not implemented")
    end()

    return ir
