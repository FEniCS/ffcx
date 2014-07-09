
from six.moves import xrange as range
from ufl.classes import Terminal

from uflacs.utils.log import error, uflacs_assert
from uflacs.analysis.datastructures import int_array, object_array, CRS, rows_to_crs
from uflacs.analysis.modified_terminals import terminal_modifier_types

def compute_dependencies(e2i, V, ignore_terminal_modifiers=True):
    if ignore_terminal_modifiers:
        terminalish = (Terminal,) + terminal_modifier_types
    else:
        terminalish = (Terminal,)

    num_rows = len(V)
    dependencies = object_array(num_rows)
    num_nonzeros = 0
    for i, v in enumerate(V):
        if isinstance(v, terminalish):
            dependencies[i] = ()
        else:
            dependencies[i] = tuple(e2i[o] for o in v.operands())
            num_nonzeros += len(dependencies[i])

    return rows_to_crs(dependencies, num_rows, num_nonzeros, int)


def mark_active(dependencies, targets):
    """Return an array marking the recursive dependencies of targets.

    Input:
    - dependencies - CRS of ints, a mapping from a symbol to the symbols of its dependencies.
    - targets      - Sequence of symbols to mark the dependencies of.

    Output:
    - active   - Truth value for each symbol.
    - num_used - Number of true values in active array.
    """
    n = len(dependencies)

    # Initial state where nothing is marked as used
    active = int_array(n)
    num_used = 0

    # Seed with initially used symbols
    active[targets] = 1

    # Mark dependencies by looping backwards through symbols array
    for s in range(n-1, -1, -1):
        if active[s]:
            num_used += 1
            active[dependencies[s]] = 1

    # Return array marking which symbols are used and the number of positives
    return active, num_used


def mark_image(inverse_dependencies, sources):
    """Return an array marking the set of symbols dependent on the sources.

    Input:
    - dependencies - CRS of ints, a mapping from a symbol to the symbols of its dependencies.
    - sources      - Sequence of symbols to mark the dependants of.

    Output:
    - image    - Truth value for each symbol.
    - num_used - Number of true values in active array.
    """
    n = len(inverse_dependencies)

    # Initial state where nothing is marked as used
    image = int_array(n)
    num_used = 0

    # Seed with initially used symbols
    image[sources] = 1

    # Mark dependencies by looping forwards through symbols array
    for s in range(n):
        if image[s]:
            num_used += 1
            image[inverse_dependencies[s]] = 1

    # Return array marking which symbols are used and the number of positives
    return image, num_used
