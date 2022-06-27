# Copyright (C) 2022 Igor A. Baratta
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import collections

from ffcx.codegeneration.indices import MultiIndex


def fuse_loops(lang, definitions):
    """
    Merge a sequence of loops with the same iteration space into a single loop.

    Loop fusion improves data locality, cache reuse and decreases the loop control overhead.
    NOTE: Loop fusion might increase the pressure on register allocation.
    Ideally, we should define a cost function to determine how many loops should fuse at a time.
    """

    bodies = collections.defaultdict(list)
    indices = collections.defaultdict(MultiIndex)

    pre_loop = []

    for access, definition in definitions.items():
        for d in definition:
            if isinstance(d, lang.NestedForRange):
                index = d.multi_indices[0]
                if index.dim == 1:
                    hash_ = hash(index)
                    bodies[hash_] += [d.body()]
                    indices[hash_] = index
                else:
                    pre_loop += [d]
            else:
                pre_loop += [d]

    fused = []
    for key in indices.keys():
        body = bodies[key]
        index = indices[key]
        fused += [lang.NestedForRange([index], body)]

    code = []
    code += pre_loop
    code += fused
    return code


def extract_indices(lang, expression):
    indices = []
    if (isinstance(expression, tuple)):
        for expr in expression:
            indices += extract_indices(lang, expr)
    elif hasattr(expression, 'args'):
        for expr in expression.args:
            indices += extract_indices(lang, expr)
    elif hasattr(expression, 'rhs'):
        indices += extract_indices(lang, expression.rhs)
        indices += extract_indices(lang, expression.lhs)
    elif isinstance(expression, lang.Symbol):
        indices += [expression]
    else:
        indices += []
    return indices


def compute_sizes(index_set, indices, sizes):
    positions = [indices.index(idx) for idx in index_set]
    ranges = [sizes[p] for p in positions]
    return ranges


# Let A, B, C be nd-tensors, this function generates code
# for computing the tensor contraction
# C{Ic} = A{Ia} * B {Ib}
def tensor_contraction(lang, A, B, C, indices, sizes, scalar_type):
    code = []

    Ib = set(extract_indices(lang, B.indices))
    Ia = set(extract_indices(lang, A.indices))

    # contracted index
    Ik = Ib.intersection(Ia)
    # Loop indices
    Iu = Ib.union(Ia)
    # output indices
    Ic = Iu.difference(Ik)

    c_sizes = compute_sizes(Ic, indices, sizes)
    code += [lang.ArrayDecl(scalar_type, C, c_sizes, values=0.0)]

    lhs = C[list(Ic)]
    rhs = lang.Mul(A, B)

    body = lang.AssignAdd(lhs, rhs)

    Iu = list(Iu)
    u_sizes = compute_sizes(Iu, indices, sizes)
    for i in reversed(range(len(Iu))):
        body = body
        for_range = lang.ForRange(Iu[i], 0, u_sizes[i], body=body)
        body = for_range

    code += [body]

    return code, lhs


def assign_add(lang, A, B, indices, sizes):
    code = []
    Ia = set(extract_indices(lang, A.indices))
    Ib = set(extract_indices(lang, B.indices))

    assert Ia == Ib, "Only copy if index sets are equal"
    Ia = list(Ia)
    a_sizes = compute_sizes(Ia, indices, sizes)

    body = lang.AssignAdd(B, A)
    for i in reversed(range(len(Ia))):
        body = body
        for_range = lang.ForRange(Ia[i], 0, a_sizes[i], body=body)
        body = for_range

    code += [body]
    return code


def sum_factorise(lang, expression, scalar_type):
    assert isinstance(expression, lang.NestedForRange)
    counter = 0
    indices = expression.indices
    sizes = expression.ranges
    if len(indices) < 4 or len(indices) > 6:
        return expression
    else:
        code = []
        expr = expression.body().expr
        assert isinstance(expr.rhs, (lang.Sum, lang.Product))

        if isinstance(expr.rhs, lang.Product):
            terms = [expr.rhs]
        else:
            terms = expr.rhs.args

        for term in terms:
            B = term.args[0]
            tables = term.args[1]
            for phi in tables.args:
                C = lang.Symbol(f"temp{counter}")
                t_code, B = tensor_contraction(lang, phi, B, C, indices, sizes, scalar_type)
                code += t_code
                counter += 1
            code += assign_add(lang, B, expr.lhs, indices, sizes)
        code = lang.Scope(code)
        return code
