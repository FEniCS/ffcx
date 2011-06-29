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
#
# Last changed: 2011-06-29

from ufl import Coefficient
from ufl.algorithms.analysis import extract_arguments

from ffc.log import info, error
from ffc.compiler import compile_form

def prepare_input_arguments(forms, object_names, reserved_objects):

    print object_names
    print reserved_objects
    # Extract unknown (of None if not defined)
    u = reserved_objects.get("unknown", None)

    # Nonlinear case
    if len(forms) == 2:
        (F, M) = forms

        #Check that unknown is defined
        assert (u), "Variable 'unknown' must be defined!"
        return (F, M, u)

    elif len(forms) == 3:
        # If unknown is undefined, define discrete solution as coefficient
        # on trial element
        (a, L, M) = forms
        if u is None:
            V = extract_arguments(a)[1].element()
            u = Coefficient(V)
            object_names[id(u)] = "__discrete_primal_solution"
        return ((a, L), M, u)
    else:
        error("Foo")

def compile_with_error_control(forms, object_names, reserved_objects,
                               prefix, parameters):

    # Check input arguments
    F, M, u = prepare_input_arguments(forms, object_names, reserved_objects)

    # Generate forms to be used for the error control
    from ffc.errorcontrol.errorcontrolgenerators import UFLErrorControlGenerator
    generator = UFLErrorControlGenerator(F, M, u)
    ec_forms = generator.generate_all_error_control_forms()

    # Check that there are no conflicts between user defined and
    # generated names
    ec_names = generator.ec_names
    comment = "%s are reserved error control names." % str(ec_names.keys())
    assert not (set(object_names.values()) & set(ec_names.values())), \
               "Conflict between user defined and generated names: %s" % comment

    # Add names generated for error control to object_names
    for (name, value) in ec_names.iteritems():
        object_names[name] = value

    # Compile error control and input (pde + goal) forms as normal
    forms = generator.primal_forms()
    compile_form(ec_forms + forms, object_names, prefix, parameters)

    return 0
