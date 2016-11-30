from ffc.log import error

import traceback
import os


try:
    from coffee.plan import ASTKernel
except ImportError:
    msg = traceback.format_exc()
    def optimize_integral_ir(ir, parameters):
        error("Failed to import coffee.plan.ASTKernel needed for optimized tsfc"
              " representation; the error message was: {}{}"
              .format(os.linesep, msg))
else:
    def optimize_integral_ir(ir, parameters):
        knl = ASTKernel(ir["tsfc"])
        knl.plan_cpu(dict(optlevel='O2'))  # TODO: optlevel from parameters
        return ir  # AST was modified in-place
