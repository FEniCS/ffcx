# Copyright (C) 2022 Igor A. Baratta
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import collections
from ffcx.codegeneration.indices import MultiIndex


def flatten_list(original_list):
    flat_list = []
    if isinstance(original_list, (list, tuple)):
        for sublist in original_list:
            flat_list += flatten_list(sublist)
        return flat_list
    else:
        return [original_list]


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
    # flatten_list = list(itertools.chain(*definitions.items()))
    #

    definitions = flatten_list(list(definitions.items()))
    definitions = [d for d in definitions if not isinstance(d, str)]

    for definition in definitions:
        if isinstance(definition, lang.NestedForRange):
            if definition.depth == 1:
                index_set = tuple(definition.indices)
                hash_ = hash(index_set)
                bodies[hash_] += [definition.body()]
                indices[hash_] = definition.multi_indices[0]
            else:
                pre_loop += [definition]
        else:
            pre_loop += [definition]

    fused = []

    for key in indices.keys():
        body = bodies[key]
        index = indices[key]
        fused += [lang.NestedForRange([index], body)]

    # for access, definition in definitions.items():
    #     for defs in definition:
    #         if not isinstance(defs, list):
    #             defs = [defs]

    #         for d in defs:
    #             if isinstance(d, lang.NestedForRange):
    #                 index = d.multi_indices[0]
    #                 if index.dim >= 1:
    #                     hash_ = hash(index)
    #                     bodies[hash_] += [d.body()]
    #                     indices[hash_] = index
    #                 else:
    #                     pre_loop += [d]
    #             else:
    #                 pre_loop += [d]

    # fused = []
    # for key in indices.keys():
    #     body = bodies[key]
    #     index = indices[key]
    #     fused += [lang.NestedForRange([index], body)]

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


def extract_multi_index(lang, tensor, all_indices, all_sizes):
    indices = extract_indices(lang, tensor.indices)
    sizes = compute_sizes(indices, all_indices, all_sizes)
    return MultiIndex(lang, indices, sizes)


def transpose_tensor(lang, A, Ia, Ib, scalar_type):
    code = []
    name = A.array.name + "transp"
    A_t = lang.Symbol(name)
    code += [lang.ArrayDecl(scalar_type, A_t, [Ib.global_size()], values=0.0)]
    body = lang.AssignAdd(A_t[Ib.global_idx()], A)
    code += [lang.NestedForRange([Ib], body)]
    return code, A_t


# Let A, B, C be nd-tensors, this function generates code
# for computing the tensor contraction
# C{Ic} = A{Ia} * B {Ib}
def tensor_contraction(lang, A, B, C, indices, sizes, scalar_type):
    code = []

    Ib = extract_multi_index(lang, B, indices, sizes)
    Ia = extract_multi_index(lang, A, indices, sizes)

    Ik = Ib.intersection(Ia)
    assert Ik.dim == 1, "contract one index at a time"

    Iu = Ia.union(Ib)
    Ic = Iu.difference(Ik)

    if isinstance(C, lang.Symbol):
        code += [lang.ArrayDecl(scalar_type, C, [Ic.global_size()], values=0.0)]

    use_gemm = True
    if use_gemm:
        Jb = Ib.intersection(Ic)
        Ib_ = Ik.union(Jb)
        Iu = Ia.union(Ib_)
        Ic_ = Iu.difference(Ik)

        transp_code, newB = transpose_tensor(lang, B, Ib, Ib_, scalar_type)

        Jb = Jb.collapse("id")
        Ib_ = Ik.union(Jb)
        Iu = Ia.union(Ib_)
        Ic_ = Iu.difference(Ik)

        code += transp_code
        rhs = lang.Mul(A, newB[Ib_.global_idx()])
        if isinstance(C, lang.Symbol):
            lhs = C[Ic_.global_idx()]
        elif isinstance(C, lang.ArrayAccess):
            lhs = C.array[Ic_.global_idx()]
        body = lang.AssignAdd(lhs, rhs)
        code += [lang.NestedForRange([Iu], body)]
        lhs = C[Ic.global_idx()]
    else:
        if isinstance(C, lang.Symbol):
            lhs = C[Ic.global_idx()]
        else:
            lhs = C
        rhs = lang.Mul(A, B)
        body = lang.AssignAdd(lhs, rhs)
        code += [lang.NestedForRange([Iu], body)]

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
            for phi in tables.args[:-1]:
                C = lang.Symbol(f"temp{counter}")
                t_code, B = tensor_contraction(lang, phi, B, C, indices, sizes, scalar_type)
                counter += 1
                code += t_code
            phi = tables.args[-1]
            t_code, B = tensor_contraction(lang, phi, B, expr.lhs, indices, sizes, scalar_type)
            code += t_code
        code = lang.Scope(code)

        return code
