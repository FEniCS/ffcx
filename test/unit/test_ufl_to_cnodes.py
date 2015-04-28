
from ufl import as_ufl
from uflacs.language.ufl_to_cnodes import ufl_to_cnodes

def test_ufl_to_cnodes():
    examples = [
        (as_ufl(1.0), "1.0"),
        (as_ufl(2), "2"),
        (as_ufl(0.0), "0.0"),
        #(as_ufl(0), "0"), # UFL doesn't preserve float type for zeros...
        (as_ufl(0.0), "0.0"),
        ]
    for expr, code in examples:
        assert str(ufl_to_cnodes(expr)) == code

def understanding_generate_partition():
    which_quadrature_loop = num_points

    def access_node(i, v):
        select on type of v:
            case literal:
                vaccess = literal
            case other modified terminal:
                vaccess = backend_access(v/mt, table_ranges[i], which_quadrature_loop)

    for i, v in partition:
        #switch type(v):
        if literal:
            vaccess = ufl_to_backend_literal(v)
            vdef = None
            vexpr = None
        elif new definition:
            vdef = backend_definitions.ufl_to_cstatement(v)
            terminalcode += vdef

            vcode = backend_language.ufl_to_cexpr()

    definitions = []
    intermediates = []

    for i in partition_indices:
        v = V[i]

        based on type(v):
            # Define how to compute v and access the result of the computation:
            vdef = None or CStatement # CStatement defining new variables
            vexpr = CExpr # CExpr with expression
            vaccess = vexpr or intermediate variable access # CExpr for future reference

            cases:
                modified terminal:
                    vdef = None or CStatement
                    vexpr = CExpr
                    vaccess = vexpr or backend defined variable access
                operator:
                    vdef = None or CStatement
                    vexpr = CExpr
                    vaccess = vexpr or intermediate variable access

        # Rename variables -> access

        # Add definition of new variable if applicable
        if vdef is not None:
            definitions.append(vdef)

        # Add computation of new intermediate if applicable
        if vaccess is not vexpr:
            intermediates.append(Assign(vaccess, vexpr))

        # Add to accesses mapping for future reference: (TODO: Which one?)
        vaccesses[v] = vaccess
        iaccesses[i] = vaccess
