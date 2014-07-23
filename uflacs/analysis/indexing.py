

# FIXME: Clean up duplicates in this module, cover with tests and profile


from six.moves import xrange as range

from ufl.common import product, sorted_by_count
from ufl.permutation import compute_indices
from ufl.utils.indexflattening import shape_to_strides, flatten_multiindex
from ufl.classes import ComponentTensor
from ufl.classes import FixedIndex
from ufl.classes import Index
from ufl.classes import Indexed


def map_indexed_arg_components(indexed):  # FIXME: This is the one in use. Is it the best?
    assert isinstance(indexed, Indexed)
    e1 = indexed
    e2, mi = e1.operands()
    d = _map_indexed_components(e2, e1, mi)
    assert all(isinstance(x, int) for x in d)
    assert len(set(d)) == len(d)
    return d


def _map_indexed_components(tensor, indexed, multiindex):
    e2 = tensor
    e1 = indexed  # e1 = e2[multiindex]

    # Compute index shape
    ind1 = sorted_by_count(e1.free_indices())
    ind2 = sorted_by_count(e2.free_indices())
    idims1 = e1.index_dimensions()
    idims2 = e2.index_dimensions()
    ish1 = tuple(idims1[i] for i in ind1)
    ish2 = tuple(idims2[i] for i in ind2)

    # Compute regular and total shape
    sh1 = e1.shape()
    sh2 = e2.shape()
    tsh1 = sh1 + ish1
    tsh2 = sh2 + ish2
    # r1 = len(tsh1)
    r2 = len(tsh2)
    # str1 = shape_to_strides(tsh1)
    str2 = shape_to_strides(tsh2)
    assert not sh1
    assert sh2  # Must have shape to be indexed in the first place
    assert product(tsh1) <= product(tsh2)

    # Build map from ind2/ish2 position (-offset nmui) to ind1/ish1 position
    ind2_to_ind1_map = [None] * len(ind2)
    for k, i in enumerate(ind2):
        ind2_to_ind1_map[k] = ind1.index(i)

    # Build map from ind1/ish1 position to mi position
    nmui = len(multiindex)
    multiindex_to_ind1_map = [None] * nmui
    for k, i in enumerate(multiindex):
        if isinstance(i, Index):
            multiindex_to_ind1_map[k] = ind1.index(i)

    # Build map from flattened e1 component to flattened e2 component
    perm1 = compute_indices(tsh1)
    ni1 = product(tsh1)

    # Situation: e1 = e2[mi]
    d1 = [None] * ni1
    p2 = [None] * r2
    assert len(sh2) == nmui
    p2ks = set()
    for k, i in enumerate(multiindex):
        if isinstance(i, FixedIndex):
            p2[k] = int(i)
            p2ks.add(k)
    for c1, p1 in enumerate(perm1):
        for k, i in enumerate(multiindex):
            if isinstance(i, Index):
                p2[k] = p1[multiindex_to_ind1_map[k]]
                p2ks.add(k)
        for k, i in enumerate(ind2_to_ind1_map):
            p2[nmui + k] = p1[i]
            p2ks.add(nmui + k)
        c2 = flatten_multiindex(p2, str2)
        d1[c1] = c2

    return d1


def map_component_tensor_arg_components(component_tensor):  # FIXME: This is the one in use. Is it the best?
    assert isinstance(component_tensor, ComponentTensor)
    e2 = component_tensor
    e1, mi = e2.operands()
    d = _map_component_tensor_components(e2, e1, mi)
    assert all(isinstance(x, int) for x in d)
    assert len(set(d)) == len(d)
    return d


def _map_component_tensor_components(tensor, indexed, multiindex):
    e1 = indexed
    e2 = tensor  # e2 = as_tensor(e1, multiindex)
    mi = [i for i in multiindex if isinstance(i, Index)]

    # Compute index shape
    ind1 = sorted_by_count(e1.free_indices())
    ind2 = sorted_by_count(e2.free_indices())
    idims1 = e1.index_dimensions()
    idims2 = e2.index_dimensions()
    ish1 = tuple(idims1[i] for i in ind1)
    ish2 = tuple(idims2[i] for i in ind2)

    # Compute regular and total shape
    sh1 = e1.shape()
    sh2 = e2.shape()
    tsh1 = sh1 + ish1
    tsh2 = sh2 + ish2
    r1 = len(tsh1)
    r2 = len(tsh2)
    str1 = shape_to_strides(tsh1)
    # str2 = shape_to_strides(tsh2)
    assert not sh1
    assert sh2
    assert len(mi) == len(multiindex)
    assert product(tsh1) == product(tsh2)
    assert ish1

    assert all(i in ind1 for i in ind2)

    nmui = len(multiindex)
    assert nmui == len(sh2)

    # Build map from ind2/ish2 position (-offset nmui) to ind1/ish1 position
    p2_to_p1_map = [None] * r2
    for k, i in enumerate(ind2):
        p2_to_p1_map[k + nmui] = ind1.index(i)

    # Build map from ind1/ish1 position to mi position
    for k, i in enumerate(mi):
        p2_to_p1_map[k] = ind1.index(mi[k])

    # Build map from flattened e1 component to flattened e2 component
    perm2 = compute_indices(tsh2)
    ni2 = product(tsh2)

    # Situation: e2 = as_tensor(e1, mi)
    d2 = [None] * ni2
    p1 = [None] * r1
    for c2, p2 in enumerate(perm2):
        for k2, k1 in enumerate(p2_to_p1_map):
            p1[k1] = p2[k2]
        c1 = flatten_multiindex(p1, str1)
        d2[c2] = c1

    return d2


def __map_indexed_to_arg_components(indexed):
    e1 = indexed
    assert isinstance(e1, Indexed)
    A1, mi1 = e1.operands()
    e2 = A1

    # Compute index shape
    ind1 = sorted_by_count(e1.free_indices())
    ind2 = sorted_by_count(e2.free_indices())
    idims1 = e1.index_dimensions()
    idims2 = e2.index_dimensions()
    ish1 = tuple(idims1[i] for i in ind1)
    ish2 = tuple(idims2[i] for i in ind2)

    # Compute regular and total shape
    sh1 = e1.shape()
    sh2 = e2.shape()
    tsh1 = sh1 + ish1
    tsh2 = sh2 + ish2
    str1 = shape_to_strides(tsh1)
    str2 = shape_to_strides(tsh2)
    assert product(tsh1) == product(tsh2)
    assert (not sh1) and (ish1) and (sh2) and (not ish2)

    sh_to_ind_map = [ind1.index(i) for i in mi1 if isinstance(i, Index)]
    comp1 = []
    comp2 = []
    for p2 in compute_indices(sh2):
        p1 = [None] * len(p2)
        for j, p in enumerate(p2):
            p1[sh_to_ind_map[j]] = p
        c1 = flatten_multiindex(p1, str1)
        c2 = flatten_multiindex(p2, str2)
        # print c1, c2
        comp1.append(c1)
        comp2.append(c2)
    return tuple(comp1), tuple(comp2)


def __map_indexed_arg_components4(indexed):
    assert isinstance(indexed, Indexed)
    e1 = indexed
    e2, mi = e1.operands()

    # Compute index shape
    ind1 = sorted_by_count(e1.free_indices())
    ind2 = sorted_by_count(e2.free_indices())
    idims1 = e1.index_dimensions()
    idims2 = e2.index_dimensions()
    ish1 = tuple(idims1[i] for i in ind1)
    ish2 = tuple(idims2[i] for i in ind2)

    # Compute regular and total shape
    sh1 = e1.shape()
    sh2 = e2.shape()
    tsh1 = sh1 + ish1
    tsh2 = sh2 + ish2
    str1 = shape_to_strides(tsh1)
    # str2 = shape_to_strides(tsh2)
    assert product(tsh1) == product(tsh2)
    assert (not sh1) and (ish1) and (sh2) and (not ish2)

    # Build map from ind1/ish1 position to mi position
    mi = [i for i in mi if isinstance(i, Index)]
    nmi = len(mi)
    ind1_to_mi_map = [None] * nmi
    for k in range(nmi):
        ind1_to_mi_map[ind1.index(mi[k])] = k

    # Build map from flattened e1 component to flattened e2 component
    indices2 = compute_indices(sh2)
    ni = len(indices2)
    d1 = [None] * ni
    d2 = [None] * ni
    for c2, p2 in enumerate(indices2):
        p1 = [p2[k] for k in ind1_to_mi_map]
        c1 = flatten_multiindex(p1, str1)
        d1[c1] = c2
        d2[c2] = c1
    assert d1 == d2
    return d1


def __map_component_tensor_arg_components4(component_tensor):
    assert isinstance(component_tensor, ComponentTensor)
    e2 = component_tensor
    e1, mi = e2.operands()

    # Compute index shape
    ind1 = sorted_by_count(e1.free_indices())
    ind2 = sorted_by_count(e2.free_indices())
    idims1 = e1.index_dimensions()
    idims2 = e2.index_dimensions()
    ish1 = tuple(idims1[i] for i in ind1)
    ish2 = tuple(idims2[i] for i in ind2)

    # Compute regular and total shape
    sh1 = e1.shape()
    sh2 = e2.shape()
    tsh1 = sh1 + ish1
    tsh2 = sh2 + ish2
    str1 = shape_to_strides(tsh1)
    # str2 = shape_to_strides(tsh2)
    assert product(tsh1) == product(tsh2)
    assert (not sh1) and (ish1) and (sh2) and (not ish2)

    # Build map from ind1/ish1 position to mi position
    mi = [i for i in mi if isinstance(i, Index)]
    nmi = len(mi)
    ind1_to_mi_map = [None] * nmi
    for k in range(nmi):
        ind1_to_mi_map[ind1.index(mi[k])] = k

    # Build map from flattened e1 component to flattened e2 component
    indices2 = compute_indices(sh2)
    ni = len(indices2)
    d1 = [None] * ni
    d2 = [None] * ni
    for c2, p2 in enumerate(indices2):
        p1 = [p2[k] for k in ind1_to_mi_map]
        c1 = flatten_multiindex(p1, str1)
        d1[c1] = c2
        d2[c2] = c1
    assert d1 == d2
    return d2
