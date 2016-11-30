import coffee
from coffee.plan import ASTKernel


def optimize_integral_ir(ir, parameters):
    knl = ASTKernel(ir["tsfc"])
    knl.plan_cpu(dict(optlevel='O2'))  # TODO: optlevel from parameters
    return ir  # AST was modified in-place
