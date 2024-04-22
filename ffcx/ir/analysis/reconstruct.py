# Copyright (C) 2011-2017 Martin Sandve Alnæs and Chris Richardson
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Reconstruct."""

import ufl


def handle_scalar_nary(o, ops):
    """Handle a scalary nary operator."""
    if o.ufl_shape != ():
        raise RuntimeError("Expecting scalar.")
    sops = [op[0] for op in ops]
    return [o._ufl_expr_reconstruct_(*sops)]


def handle_condition(o, ops):
    """Handle a condition."""
    # A condition is always scalar, so len(op) == 1
    sops = [op[0] for op in ops]
    return [o._ufl_expr_reconstruct_(*sops)]


def handle_conditional(o, ops):
    """Handle a conditional."""
    # A condition can be non scalar
    symbols = []
    n = len(ops[1])
    if len(ops[0]) != 1:
        raise RuntimeError("Condition should be scalar.")
    if n != len(ops[2]):
        raise RuntimeError("Conditional branches should have same shape.")
    for i in range(len(ops[1])):
        sops = (ops[0][0], ops[1][i], ops[2][i])
        symbols.append(o._ufl_expr_reconstruct_(*sops))
    return symbols


def handle_elementwise_unary(o, ops):
    """Handle a elementwise unary operator."""
    if len(ops) > 1:
        raise RuntimeError("Expecting unary operator.")
    return [o._ufl_expr_reconstruct_(op) for op in ops[0]]


def handle_division(o, ops):
    """Handle a division."""
    if len(ops) != 2:
        raise RuntimeError("Expecting two operands.")
    if len(ops[1]) != 1:
        raise RuntimeError("Expecting scalar divisor.")
    (b,) = ops[1]
    return [o._ufl_expr_reconstruct_(a, b) for a in ops[0]]


def handle_sum(o, ops):
    """Handle a sum."""
    if len(ops) != 2:
        raise RuntimeError("Expecting two operands.")
    if len(ops[0]) != len(ops[1]):
        raise RuntimeError("Expecting scalar divisor.")
    return [o._ufl_expr_reconstruct_(a, b) for a, b in zip(ops[0], ops[1])]


def handle_product(o, ops):
    """Handle a product."""
    if len(ops) != 2:
        raise RuntimeError("Expecting two operands.")

    # Get the simple cases out of the way
    if len(ops[0]) == 1:  # True scalar * something
        (a,) = ops[0]
        return [ufl.classes.Product(a, b) for b in ops[1]]
    elif len(ops[1]) == 1:  # Something * true scalar
        (b,) = ops[1]
        return [ufl.classes.Product(a, b) for a in ops[0]]

    # Neither of operands are true scalars, this is the tricky part
    o0, o1 = o.ufl_operands

    # Get shapes and index shapes
    fi = o.ufl_free_indices
    fi0 = o0.ufl_free_indices
    fi1 = o1.ufl_free_indices
    fid = o.ufl_index_dimensions
    fid0 = o0.ufl_index_dimensions
    fid1 = o1.ufl_index_dimensions

    # Need to map each return component to one component of o0 and
    # one component of o1
    indices = ufl.permutation.compute_indices(fid)

    # Compute which component of o0 is used in component (comp,ind) of o
    # Compute strides within free index spaces
    ist0 = ufl.utils.indexflattening.shape_to_strides(fid0)
    ist1 = ufl.utils.indexflattening.shape_to_strides(fid1)
    # Map o0 and o1 indices to o indices
    indmap0 = [fi.index(i) for i in fi0]
    indmap1 = [fi.index(i) for i in fi1]
    indks = [
        (
            ufl.utils.indexflattening.flatten_multiindex([ind[i] for i in indmap0], ist0),
            ufl.utils.indexflattening.flatten_multiindex([ind[i] for i in indmap1], ist1),
        )
        for ind in indices
    ]

    # Build products for scalar components
    results = [ufl.classes.Product(ops[0][k0], ops[1][k1]) for k0, k1 in indks]
    return results


def handle_index_sum(o, ops):
    """Handle an index sum."""
    summand, mi = o.ufl_operands
    ic = mi[0].count()
    fi = summand.ufl_free_indices
    fid = summand.ufl_index_dimensions
    ipos = fi.index(ic)
    d = fid[ipos]

    # Compute "macro-dimensions" before and after i in the total shape of a
    predim = ufl.product(summand.ufl_shape) * ufl.product(fid[:ipos])
    postdim = ufl.product(fid[ipos + 1 :])

    # Map each flattened total component of summand to
    # flattened total component of indexsum o by removing
    # axis corresponding to summation index ii.
    ss = ops[0]  # Scalar subexpressions of summand
    if len(ss) != predim * postdim * d:
        raise RuntimeError("Mismatching number of subexpressions.")
    sops = []
    for i in range(predim):
        iind = i * (postdim * d)
        for k in range(postdim):
            ind = iind + k
            sops.append([ss[ind + j * postdim] for j in range(d)])

    # For each scalar output component, sum over collected subcomponents
    # TODO: Need to split this into binary additions to work with future CRSArray format,
    #       i.e. emitting more expressions than there are symbols for this node.
    results = [sum(sop) for sop in sops]
    return results


# TODO: To implement compound tensor operators such as dot and inner,
# we need to identify which index to do the contractions over,
# and build expressions such as sum(a*b for a,b in zip(aops, bops))


_reconstruct_call_lookup = {
    ufl.classes.MathFunction: handle_scalar_nary,
    ufl.classes.Abs: handle_scalar_nary,
    ufl.classes.MinValue: handle_scalar_nary,
    ufl.classes.MaxValue: handle_scalar_nary,
    ufl.classes.Real: handle_elementwise_unary,
    ufl.classes.Imag: handle_elementwise_unary,
    ufl.classes.Power: handle_scalar_nary,
    ufl.classes.BesselFunction: handle_scalar_nary,
    ufl.classes.Atan2: handle_scalar_nary,
    ufl.classes.Product: handle_product,
    ufl.classes.Division: handle_division,
    ufl.classes.Sum: handle_sum,
    ufl.classes.IndexSum: handle_index_sum,
    ufl.classes.Conj: handle_elementwise_unary,
    ufl.classes.Conditional: handle_conditional,
    ufl.classes.Condition: handle_condition,
}


def reconstruct(o, *args):
    """Reconstruct."""
    # First look for exact match
    f = _reconstruct_call_lookup.get(type(o), False)
    if f:
        return f(o, *args)
    else:
        # Look for parent class types instead
        for k in _reconstruct_call_lookup.keys():
            if isinstance(o, k):
                return _reconstruct_call_lookup[k](o, *args)
        # Nothing found
        raise RuntimeError("Not expecting expression of type %s in here." % type(o))
