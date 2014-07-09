
import numpy as np

default_tolerance = 1e-14

def equal_tables(a, b, eps=default_tolerance):
    "Compare tables to be equal within a tolerance."
    a = np.asarray(a)
    b = np.asarray(b)
    if a.shape != b.shape:
        return False
    if len(a.shape) > 1:
        return all(equal_tables(a[i], b[i], eps) for i in xrange(a.shape[0]))
    def scalars_equal(x, y, eps):
        return abs(x-y) < eps
    return all(scalars_equal(a[i], b[i], eps) for i in xrange(a.shape[0]))

def strip_table_zeros(table, eps=default_tolerance):
    "Strip zero columns from table. Returns column range (begin,end) and the new compact table."
    # Get shape of table and number of columns, defined as the last axis
    table = np.asarray(table)
    sh = table.shape
    nc = sh[-1]

    # Find first nonzero column
    begin = nc
    for i in xrange(nc):
        if np.linalg.norm(table[...,i]) > eps:
            begin = i
            break

    # Find (one beyond) last nonzero column
    end = begin
    for i in xrange(nc-1,begin-1,-1):
        if np.linalg.norm(table[...,i]) > eps:
            end = i+1
            break

    # Make subtable by stripping first and last columns
    return begin, end, table[...,begin:end]

def build_unique_tables(tables):
    """Given a list or dict of tables, return a list of unique tables
    and a dict of unique table indices for each input table key."""
    unique = []
    mapping = {}
    if isinstance(tables, list):
        keys = list(range(len(tables)))
    elif isinstance(tables, dict):
        keys = sorted(tables.keys())
    for k in keys:
        t = tables[k]
        found = -1
        for i,u in enumerate(unique):
            if equal_tables(u, t):
                found = i
                break
        if found == -1:
            i = len(unique)
            unique.append(t)
        mapping[k] = i
    return unique, mapping

def build_unique_tables2(tables):
    """Given a list or dict of tables, return a list of unique tables
    and a dict of unique table indices for each input table key."""
    unique = []
    mapping = {}

    if isinstance(tables, list):
        keys = list(range(len(tables)))
    elif isinstance(tables, dict):
        keys = sorted(tables.keys())

    for k in keys:
        t = tables[k]
        found = -1
        for i,u in enumerate(unique):
            if equal_tables(u, t):
                found = i
                break
        if found == -1:
            i = len(unique)
            unique.append(t)
        mapping[k] = i

    return unique, mapping

def get_ffc_table_values(tables, entitytype, num_points, element, flat_component, derivative_counts):
    """Extract values from ffc element table.

    Returns a 3D numpy array with axes
    (entity number, quadrature point number, dof number)
    """
    # Get quadrule/element subtable
    element_table = tables[num_points][element]

    # Temporary fix for new table structure TODO: Handle avg properly
    if len(element_table) != 1:
        print; print element_table
    assert len(element_table) == 1
    element_table = element_table[None]

    # FFC property:
    #element_counter = element_map[num_points][element]

    # Figure out shape of final array by inspecting tables
    num_entities = len(element_table)
    num_dofs = len(element_table.itervalues().next()[derivative_counts])

    # Make 3D array for final result
    shape = (num_entities, num_points, num_dofs)
    res = np.zeros(shape)

    # Loop over entities and fill table blockwise (each block = points x dofs)
    sh = element.value_shape()
    for entity in xrange(num_entities):
        # Access subtable
        entity_key = None if entitytype == "cell" else entity
        tbl = element_table[entity_key][derivative_counts]

        # Extract array for right component and order axes as (points, dofs)
        if sh == ():
            arr = np.transpose(tbl)
        else:
            arr = np.transpose(tbl[:,flat_component,:])

        # Assign block of values for this entity
        res[entity, :, :] = arr
    return res

def generate_psi_table_name(element_counter, flat_component, derivative_counts, averaged, entitytype):
    """Generate a name for the psi table of the form:
    FE#_C#_D###[_AC|_AF|][_F|V], where '#' will be an integer value.

    FE  - is a simple counter to distinguish the various bases, it will be
          assigned in an arbitrary fashion.

    C   - is the component number if any (this does not yet take into account
          tensor valued functions)

    D   - is the number of derivatives in each spatial direction if any.
          If the element is defined in 3D, then D012 means d^3(*)/dydz^2.

    AC  - marks that the element values are averaged over the cell

    AF  - marks that the element values are averaged over the facet

    F   - marks that the first array dimension enumerates facets on the cell

    v   - marks that the first array dimension enumerates vertices on the cell

    """

    name = "FE%d" % element_counter

    if isinstance(flat_component, int):
        name += "_C%d" % flat_component
    else:
        assert flat_component is None

    if any(derivative_counts):
        name += "_D" + "".join(map(str,derivative_counts))

    if averaged == "cell":
        name += "_AC"
    elif averaged == "facet":
        name += "_AF"

    if entitytype == "cell":
        pass
    elif entitytype == "facet":
        name += "_F"
    elif entitytype == "vertex":
        name += "_V"
    else:
        error("Unknown entity type %s." % entitytype)

    return name


def derivative_counts_to_listing(derivative_counts, gdim):
    derivatives = []
    for i, d in enumerate(derivative_counts):
        derivatives.extend((i,)*d)
    return derivatives

def derivative_listing_to_counts(derivatives, gdim):
    derivative_counts = [0]*gdim
    for d in derivatives:
        derivative_counts[d] += 1
    return derivative_counts

#from ufl.common import component_to_index
from ufl.permutation import build_component_numbering
def flatten_component(component, shape, symmetry):
    if shape:
        # Map component to flat index
        vi2si, si2vi = build_component_numbering(shape, symmetry)
        return vi2si[component]
    else:
        return None


def _examples(tables):
    name = generate_psi_table_name(counter, flat_component, derivative_counts, averaged, entitytype)
    values = get_ffc_table_values(tables, entitytype, num_points, element, flat_component, derivative_counts)

    begin, end, table = strip_table_zeros(table)
    all_zeros = table.shape[-1] == 0
    all_ones = equal_tables(table, np.ones(table.shape))
