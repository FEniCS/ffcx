
from six.moves import zip
from six.moves import xrange as range
from ufl.common import product
from ufl.permutation import compute_indices

import ufl
from ufl import as_vector
from ufl.classes import (MultiIndex, ComponentTensor, ListTensor, Transposed, Variable,
                         IndexSum, UtilityType, Label, ExprList, ExprMapping)
from ufl.algorithms import MultiFunction
from ufl.utils.indexflattening import flatten_multiindex, shape_to_strides
from ufl.utils.sorting import sorted_by_count

from ffc.log import error, ffc_assert
from uflacs.datastructures.arrays import int_array, object_array
from uflacs.datastructures.crs import CRS, rows_to_crs, rows_dict_to_crs
from uflacs.analysis.modified_terminals import is_modified_terminal

class ReconstructScalarSubexpressions(MultiFunction):
    def __init__(self):
        super(ReconstructScalarSubexpressions, self).__init__()

    def expr(self, o, *args, **kwargs):
        error("No handler for type %s" % type(o))

    def terminal(self, o):
        error("Not expecting terminal expression in here, got %s." % type(o))

    def scalar_nary(self, o, ops):
        ffc_assert(o.shape() == (), "Expecting scalar.")
        sops = [op[0] for op in ops]
        return [o.reconstruct(*sops)]

    # Unary scalar functions
    math_function = scalar_nary
    abs = scalar_nary
    # Binary scalar functions
    power = scalar_nary
    bessel_function = scalar_nary # TODO: Is this ok?

    def condition(self, o, ops):
        sops = [op[0] for op in ops]
        return [o.reconstruct(*sops)]

    def conditional(self, o, ops):
        sops = [op[0] for op in ops]
        return [o.reconstruct(*sops)]

    def division(self, o, ops):
        ffc_assert(len(ops) == 2, "Expecting two operands.")
        ffc_assert(len(ops[1]) == 1, "Expecting scalar divisor.")
        b, = ops[1]
        return [o.reconstruct(a, b) for a in ops[0]]

    def sum(self, o, ops):
        ffc_assert(len(ops) == 2, "Expecting two operands.")
        ffc_assert(len(ops[0]) == len(ops[1]), "Expecting scalar divisor.")
        return [o.reconstruct(a, b) for a, b in zip(ops[0], ops[1])]

    def product(self, o, ops):
        ffc_assert(len(ops) == 2, "Expecting two operands.")

        # Get the simple cases out of the way
        na = len(ops[0])
        nb = len(ops[1])

        if na == 1: # True scalar * something
            a, = ops[0]
            return [o.reconstruct(a, b) for b in ops[1]]

        if nb == 1: # Something * true scalar
            b, = ops[1]
            return [o.reconstruct(a, b) for a in ops[0]]

        # Neither of operands are true scalars, this is the tricky part
        o0, o1 = o.operands()

        def _compute_shapes(expr):
            fi = sorted_by_count(expr.free_indices())
            idims = expr.index_dimensions()
            sh = expr.shape()
            ish = tuple(idims[i] for i in fi)
            return fi, idims, sh, ish

        # Get shapes and index shapes
        fi, idims, sh, ish = _compute_shapes(o)
        fi0, idims0, sh0, ish0 = _compute_shapes(o0)
        fi1, idims1, sh1, ish1 = _compute_shapes(o1)

        # Need to map each return component to one component of o0 and one component of o1.
        components = compute_indices(sh)
        indices = compute_indices(ish)

        # Map component comp of o to component offset of o0 and o1
        if sh0:
            # o0 has shape
            im0 = product(ish0)
            st0 = shape_to_strides(sh0)
            compks = [(im0*flatten_multiindex(comp, st0), 0) for comp in components]
        elif sh1:
            # o1 has shape
            im1 = product(ish1)
            st1 = shape_to_strides(sh1)
            compks = [(0, im1*flatten_multiindex(comp, st1)) for comp in components]
        else:
            # Neither has shape (indices only)
            # (It never happens that both have shape)
            compks = [(0, 0)]

        # Compute which component of o0 is used in component (comp,ind) of o
        if fi0 or fi1:
            # Compute strides within free index spaces
            ist0 = shape_to_strides(ish0)
            ist1 = shape_to_strides(ish1)
            # Map o0 and o1 indices to o indices
            indmap0 = [fi.index(i) for i in fi0]
            indmap1 = [fi.index(i) for i in fi1]
            indks = [(flatten_multiindex([ind[i] for i in indmap0], ist0),
                      flatten_multiindex([ind[i] for i in indmap1], ist1))
                    for ind in indices]
        else:
            indks = [(0,0)]*len(indices)

        #from IPython.core.debugger import Pdb; pdb = Pdb(); pdb.set_trace()

        # Build products for scalar components
        results = []
        for k00, k10 in compks:
            for k01, k11 in indks:
                results.append(o.reconstruct(ops[0][k00 + k01], ops[1][k10 + k11]))

        #results = [o.reconstruct(ops[0][k00 + k01], ops[1][k10 + k11])
        #           for k00, k10 in compks
        #           for k01, k11 in indks]

        return results

    def index_sum(self, o, ops):
        summand, mi = o.operands()
        ii = mi[0]
        fi = summand.free_indices()
        idims = summand.index_dimensions()
        d = idims[ii]

        # Compute "macro-dimensions" before and after i in the total shape of a
        predim = product(summand.shape())
        postdim = 1
        ic = ii.count()
        for jj in fi:
            jc = jj.count()
            if jc < ic:
                predim *= idims[jj]
            elif jc > ic:
                postdim *= idims[jj]

        # Map each flattened total component of summand to
        # flattened total component of indexsum o by removing
        # axis corresponding to summation index ii.
        ss = ops[0] # Scalar subexpressions of summand
        ffc_assert(len(ss) == predim*postdim*d, "Mismatching number of subexpressions.")
        sops = []
        for i in range(predim):
            iind = i*(postdim*d)
            for k in range(postdim):
                ind = iind + k
                sops.append([ss[ind + j*postdim] for j in range(d)])

        # For each scalar output component, sum over collected subcomponents
        return [sum(sop) for sop in sops]

    # TODO: To implement compound tensor operators such as dot and inner,
    # we need to identify which index to do the contractions over,
    # and build expressions such as sum(a*b for a,b in zip(aops, bops))

def rebuild_scalar_e2i(G, DEBUG=False):
    """Build a new expression2index mapping where each subexpression is scalar valued.

    Input:
    - G.e2i
    - G.V
    - G.V_symbols
    - G.total_unique_symbols

    Output:
    - NV   - Array with reverse mapping from index to expression
    - nvs  - Tuple of ne2i indices corresponding to the last vertex of G.V

    Old output now no longer returned but possible to restore if needed:
    - ne2i - Mapping from scalar subexpressions to a contiguous unique index
    - W    - Array with reconstructed scalar subexpressions for each original symbol
    - terminals - Set of modified terminal expressions (terminals possibly wrapped
                  in grad, restriction and indexed)
    """

    # Data structures
    ne2i = {}
    NV = object_array(G.total_unique_symbols)
    W = object_array(G.total_unique_symbols)
    terminals = set()

    # These types are not expected to be part of the graph at this point
    unexpected = (MultiIndex, ExprList, ExprMapping, UtilityType, Label, ComponentTensor, ListTensor, Transposed, Variable)

    def emit_expression(s, u):
        # Allocate count for scalar expression and
        # store in all cross referenced data structures
        j = ne2i.get(u)
        if j is None:
            j = len(ne2i)
            ne2i[u] = j
            NV[j] = u
        W[s] = u
        emit_expression.W_len = max(emit_expression.W_len, s+1)
        if DEBUG: print('emitted s, j, u:', s, j, u)
    emit_expression.W_len = 0

    reconstruct_scalar_subexpressions = ReconstructScalarSubexpressions()

    handled_symbols = int_array(G.total_unique_symbols)
    for i, v in enumerate(G.V):
        # Find symbols of v components
        vs = G.V_symbols[i]

        if DEBUG: print('\n\n:: i, v, vs ::', i, v, vs)

        # Skip if there's nothing new here (should be the case for indexing types)
        if all(handled_symbols[s] for s in vs):
            continue

        #if all(W[s] is not None for s in vs):
        #    continue

        # This takes a lot of time:
        #if isinstance(v, unexpected):
        #    error("Not expecting a %s here!" % type(v))

        for s in vs:
            handled_symbols[s] = 1

        if is_modified_terminal(v):
            if 0: print("Adding terminal: ", repr(v))
            sh = v.shape()
            ffc_assert(v.free_indices() == (), "Expecting no free indices.")
            if sh == ():
                # Store single modified terminal expression component
                ffc_assert(len(vs) == 1, "Expecting single symbol for scalar valued modified terminal.")
                s, u = vs[0], v
                emit_expression(s, u)
                terminals.add(v)
            else:
                # Store each terminal expression component
                for s, c in zip(vs, compute_indices(sh)):
                    u = v[c]
                    emit_expression(s, u)
                    # FIXME: Keep modified terminal expression components in the graph that is input here!
                    terminals.add(u)
        else:
            # Find symbols of operands
            sops = []
            for j, vop in enumerate(v.operands()):
                if isinstance(vop, MultiIndex):
                    if not isinstance(v, IndexSum):
                        error("FIXME: Not expecting a %s." % type(v))
                    so = ()
                else:
                    so = G.V_symbols[G.e2i[vop]]
                sops.append(so)

            # Fetch reconstructed operand expressions
            wops = [tuple(W[k] for k in so) for so in sops]

            # Reconstruct scalar subexpressions of v
            w = reconstruct_scalar_subexpressions(v, wops)

            # Store all scalar subexpressions for v symbols
            if not len(vs) == len(w):
                print()
                print(type(v))
                print(v.shape())
                print(v.free_indices())
                print(len(vs))
                print(len(w))
                print()
                ffc_assert(len(vs) == len(w), "Expecting one symbol for each expression.")
            for s, u in zip(vs, w):
                emit_expression(s, u)

    # Reduce size of NV to the actually used parts
    ffc_assert(all(x is None for x in NV[len(ne2i):]),
                  "Expecting last part of NV to be empty.")
    NV = NV[:len(ne2i)]

    # Find symbols of final v
    vs = G.V_symbols[G.nv-1]
    ffc_assert(all(handled_symbols[s] for s in vs),
                  "Expecting that all symbols in vs are handled at this point.")
    nvs = [ne2i[W[s]] for s in vs]

    # TODO: Make it so that expressions[k] <-> NV[nvs[k][:]], len(nvs[k]) == value_size(expressions[k])

    return NV, nvs

def rebuild_expression_from_graph(G, DEBUG=False):
    "This is currently only used by tests."
    NV, nvs = rebuild_scalar_e2i(G, DEBUG=DEBUG)

    # Find expressions of final v
    w = [NV[k] for k in nvs]
    if len(w) == 1:
        return w[0]
    else:
        return as_vector(w) # TODO: Consider shape of initial v
