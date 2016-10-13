# -*- coding: utf-8 -*-
"""
This module provides compilation of forms required for goal-oriented
error control
"""

# Copyright (C) 2010 Marie E. Rognes
#
# This file is part of FFC.
#
# FFC is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# FFC is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with FFC. If not, see <http://www.gnu.org/licenses/>.

from ufl.utils.sorting import sorted_by_key
from ufl import Coefficient

from ffc.log import error
from ffc.compiler import compile_form

__all__ = ["compile_with_error_control"]


def compile_with_error_control(forms, object_names, reserved_objects,
                               prefix, parameters):
    """
    Compile forms and additionally generate and compile forms required
    for performing goal-oriented error control

    For linear problems, the input forms should be a bilinear form (a)
    and a linear form (L) specifying the variational problem and
    additionally a linear form (M) specifying the goal functional.

    For nonlinear problems, the input should be linear form (F) and a
    functional (M) specifying the goal functional.

    *Arguments*

        forms (tuple)
            Three (linear case) or two (nonlinear case) forms
            specifying the primal problem and the goal

        object_names (dict)
            Map from object ids to object names

        reserved_names (dict)
            Map from reserved object names to object ids

        prefix (string)
            Basename of header file

        parameters (dict)
            Parameters for form compilation
    """

    # Check input arguments
    F, M, u = prepare_input_arguments(forms, object_names, reserved_objects)

    # Generate forms to be used for the error control
    from ffc.errorcontrol.errorcontrolgenerators import UFLErrorControlGenerator
    generator = UFLErrorControlGenerator(F, M, u)
    ec_forms = generator.generate_all_error_control_forms()

    # Check that there are no conflicts between user defined and
    # generated names
    ec_names = generator.ec_names
    if set(object_names.values()) & set(ec_names.values()):
        comment = "%s are reserved error control names." % str(sorted(ec_names.values()))
        error("Conflict between user defined and generated names: %s" % comment)

    # Add names generated for error control to object_names
    for (objid, name) in sorted_by_key(ec_names):
        object_names[objid] = name

    # Compile error control and input (pde + goal) forms as normal
    forms = generator.primal_forms()
    code_h, code_c = compile_form(ec_forms + forms, object_names, prefix, parameters)

    return code_h, code_c


def prepare_input_arguments(forms, object_names, reserved_objects):
    """
    Extract required input arguments to UFLErrorControlGenerator.

    *Arguments*

        forms (tuple)
            Three (linear case) or two (nonlinear case) forms
            specifying the primal problem and the goal

        object_names (dict)
            Map from object ids to object names

        reserved_names (dict)
            Map from reserved object names to object ids

    *Returns*

        tuple (of length 3) containing

        Form or tuple
            A single linear form or a tuple of a bilinear and a linear
            form

        Form
            A linear form or a functional for the goal functional

        Coefficient
            The coefficient considered as the unknown
    """

    # Check that we get a tuple of forms
    if not isinstance(forms, (list, tuple)):
        error("Expecting tuple of forms, got %s" % str(forms))

    def __is_nonlinear(forms):
        return len(forms) == 2

    def __is_linear(forms):
        return len(forms) == 3

    # Extract Coefficient labelled as 'unknown'
    u = reserved_objects.get("unknown", None)

    if __is_nonlinear(forms):
        # Check that unknown is defined
        if u is None:
            error("Can't extract 'unknown'. The Coefficient representing the unknown must be labelled by 'unknown' for nonlinear problems.")

        (F, M) = forms

        # Check that forms have the expected rank
        assert len(F.arguments()) == 1
        assert len(M.arguments()) == 0

        # Return primal, goal and unknown
        return (F, M, u)

    elif __is_linear(forms):
        # Throw error if unknown is given, don't quite know what to do
        # with this case yet
        if u is not None:
            error("'unknown' defined: not implemented for linear problems")

        (a, L, M) = forms

        # Check that forms have the expected rank
        arguments = a.arguments()
        assert len(arguments) == 2
        assert len(L.arguments()) == 1
        assert len(M.arguments()) == 1

        # Standard case: create default Coefficient in trial space and
        # label it __discrete_primal_solution
        V = arguments[1].ufl_function_space()
        u = Coefficient(V)
        object_names[id(u)] = "__discrete_primal_solution"

        return ((a, L), M, u)

    else:
        error("Wrong input tuple length: got %s, expected 2 or 3-tuple"
              % str(forms))
