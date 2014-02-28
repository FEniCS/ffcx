
class ToyCppStatementFormatter(object):
    """Class containing functions for generating definitions of registers,
    argument loops, and output variable names."""
    def __init__(self, unused_dependency_handler, ir):
        self._unused_dependency_handler = unused_dependency_handler

    def define_registers(self, num_registers, partition=None):
        # TODO: Partition is not yet used by compiler
        code = ["// Declaring variables for intermediate computations:"]
        name = "s" if partition is None else ("s%d" % partition)
        code += ["double %s[%d];" % (name, num_registers)]
        code.append("")
        return code

    def define_piecewise_geometry(self):
        code = ["// TODO: Compute piecewise constant geometry here"]
        return code

    def define_coord_loop(self):
        return None

    def define_coord_vars(self):
        code = ["// TODO: Computing coordinates in necessary coordinate systems:"]
        return code

    def define_coord_dependent_geometry(self):
        code = ["// TODO: Computing coordinate dependent geometry here"]
        return code

    def accumulation_scaling_factor(self):
        return None

    def define_piecewise_coefficients(self):
        code = ["// TODO: Compute piecewise constant coefficients here"]
        return code

    def define_coord_dependent_coefficients(self):
        return ["// Compute x dependent coefficients and evt. their derivatives"]

    def define_argument_for_loop(self, argument_count):
        iname = "i%d" % (argument_count,)
        isize = "n%d" % (argument_count,)
        return "for (int %s = 0; %s < %s; ++%s)" % (iname, iname, isize, iname)

    def define_argument_loop_vars(self, argument_count):
        return ["// Compute argument %d and evt. its derivatives" % (argument_count,)]

    def define_output_variables_reset(self):
        return []

    def output_variable_names(self, num_variables):
        return ["A[%d]" % (i,) for i in xrange(num_variables)]
