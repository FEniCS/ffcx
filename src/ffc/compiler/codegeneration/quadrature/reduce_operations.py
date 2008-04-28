#!/usr/bin/env python

def operation_count(expression, format):
    # Cheat to get character for add and multiply
    add   = format["add"](["", ""])
    mult  = format["multiply"](["", ""])
    adds = len(get_products(expression, format, add)) - 1
    return expression.count(mult) + adds

def get_geo_terms(expression, geo_terms, format):
    # Cheat to get character for add and multiply
    add   = format["add"](["", ""])
    mult  = format["multiply"](["", ""])
    access = format["array access"]("")
    group = format["grouping"]("")
    l = access[0]
    r = access[1]
    num_geo = len(geo_terms)
    prods = get_products(expression, format, add)
    new_prods = []

    G = "G"
#    print "\ngeo_terms: ", geo_terms
#    print "\nexpr: ", expression
    for p in prods:
#        print "\np: ", p
        vrs = get_products(p, format, mult)
        geos = [v for v in vrs if v and not l in v and not r in v]
        geos.sort()
        geo = mult.join(geos)
#        print "geo: ", geo
#        print "vrs: ", vrs
        if geo:
            for g in geos:
                vrs.remove(g)
#                print "g: ", g
#                print "vrs: ", vrs
            if not geo in geo_terms:
                geo_terms[geo] = G + str(num_geo)
                num_geo += 1
            vrs.append(geo_terms[geo])
        new_prods.append(mult.join(vrs))

#    print "geo_terms: ", geo_terms

    return add.join(new_prods)

def get_products(expression, format, operator):

    access = format["array access"]("")
    group = format["grouping"]("")
    la = access[0]
    ra = access[1]
    lg = group[0]
    rg = group[1]

    prods = expression.split(operator)
    new_prods = [prods[0]]
    prods.pop(0)
    while prods:
        p = prods[0]
#        print "p: ", p
        if not new_prods[-1].count(la) == new_prods[-1].count(ra):
            new_prods[-1] = operator.join([new_prods[-1], p])
        elif not new_prods[-1].count(lg) == new_prods[-1].count(rg):
            new_prods[-1] = operator.join([new_prods[-1], p])
        else:
            new_prods.append(p)
        prods.remove(p)
    if expression == operator.join(new_prods):
#        print new_prods
        return new_prods
    else:
        raise RuntimeError, "Something wrong with expression"

def remove_vars(expr, vrs, format):
    mult  = format["multiply"](["", ""])

    mults = get_products(expr, format, mult)
    if not isinstance(vrs, list):
        vrs = [vrs]
    for v in vrs:
        mults.remove(v)
#    print "Rem mults: ", mults
    if not mults:
#        mults.append(str(1))
        return str(1.0)
#    print "join mults: ", mult.join(mults)
    collect = ""
    try:
        collect = str(eval(mult.join(mults)))
    except:
        collect = collect_floats(mult.join(mults), format)
   
    return collect

def get_all_variables(expression, format):
    add   = format["add"](["", ""])
    sub   = format["subtract"](["", ""])
    mult  = format["multiply"](["", ""])
    group = format["grouping"]("")

    prods = get_products(expression, format, add)
    prods = [p for p in prods if p]
#    print "sums: ", sums
    variables = []
    for i in range(len(prods)):
        p = prods[i]
        # only extract unique variables
        vrs = get_products(p, format, mult)
        for v in vrs:
            try:
                float(v)
            except:
                variables.append(v)
    return (prods, variables)

def collect_floats(expression, format):
    add   = format["add"](["", ""])
    mult  = format["multiply"](["", ""])
    div  = format["division"]
    group = format["grouping"]("")
    l = group[0]
    r = group[1]

    if div in expression:
        return expression
    if l in expression or r in expression:
        raise RuntimeError, "Illegal expression for collect floats"

    floats = 0.0
    new_prods = []
    for p in get_products(expression, format, add):
#        print "p: ", p
        new_mults = []
        mults = [m for m in get_products(p, format, mult) if m]
        fact = 1.0
        for m in mults:
            try:
#                print "m: ", m
                fact *= float(m)
            except:
                new_mults.append(m)
        new_mults.sort()
#        print "\nnew_mults: ", new_mults
        if not fact == 1.0:
            new_mults = [str(fact)] + new_mults
        p = mult.join(new_mults)
#        print "p: ", p
        try:
            floats += float(eval(p))
        except:
            new_prods.append(p)
    if floats:
        new_prods.append(str(floats))
    return add.join(new_prods)

def group_vars(expr, format):
    add   = format["add"](["", ""])
    sub   = format["subtract"](["", ""])
    mult  = format["multiply"](["", ""])
    div  = format["division"]
    group = format["grouping"]("")
    l = group[0]
    r = group[1]
    ab   = format["absolute value"]("").split(l)[0]
    sq   = format["sqrt"]("").split(l)[0]

#    print ab, sq

    if div in expr or sq in expr or ab in expr:
        return expr
    if l in expr or r in expr:
        raise RuntimeError, "Illegal expression for group_vars"
    if not expr:
        return expr
    new_prods = {}
    prods = get_products(expr, format, add)
#    print "PRODS: ", prods
    for p in prods:
        vrs = get_all_variables(p, format)
#        print "vrs: ", vrs
        var = [v for v in vrs[1] if v]
        var.sort()
        factor = remove_vars(p, var, format)
#        print "factor: ", factor
        var = mult.join(var)
#        print "var: ", var
        if var in new_prods:
            new_prods[var] += [factor]
        else:
            new_prods[var] = [factor]
    
    prods = []
#    print "new t: ", new_prods
    for p in new_prods:
#        print "p: ", p
#        f = collect_floats(add.join(new_prods[p]), format)
        f = str(eval(add.join(new_prods[p])))
        if eval(f) == 1.0:
            f = ""
#        print "f: ", f
        if f and p:
            prods.append(mult.join([f,p]))
        elif f:
            prods.append(f)
        elif p:
            prods.append(p)

#        print "terms: ", terms
    return add.join(prods)


def get_variables(expression, format):
    add   = format["add"](["", ""])
    mult  = format["multiply"](["", ""])
    group = format["grouping"]("")

    prods = get_products(expression, format, add)
    prods = [p for p in prods if p]
#    new_prods = []
#    for p in prods:
#        new_prods += get_products(p, format, mult)

#    print "prods: ", prods
#    print "new_prods: ", new_prods

    variables = {}
    for i in range(len(prods)):
        p = prods[i]
        # only extract unique variables
#        print "get rpds: ", get_products(p, format, mult)
        vrs = list(set( get_products(p, format, mult) ))
#        vrs = list(set( p.split(mult) ))
        for v in vrs:
            if v in variables:
                variables[v][0] += 1
                variables[v][1].append(i)
            else:
                variables[v] = [1, [i]]

    return (prods, variables)

def reduction_possible(variables):
    max_val = 1
    max_var = ""
    max_vars = []
    for key in variables:
        if max_val < variables[key][0]:
            max_val = variables[key][0]
            max_var = key

    # If we found a variable that appears in products multiple times, check if
    # other variables appear in the exact same products
    if max_var:
        for key in variables:
            # Check if we have more variables in the same products
            if max_val == variables[key][0] and variables[max_var][1] == variables[key][1]:
                max_vars.append(key)
#    print "max_val: ", max_val
#    print "max_var: ", max_var
    return max_vars

def reduce_operations(expression, format):

    # Cheat to get character for add and multiply
    add   = format["add"](["", ""])
    mult  = format["multiply"](["", ""])
    group = format["grouping"]("")

#    print expression
    expression = expand_operations(expression, format)
#    print expression

    prods, variables = get_variables(expression, format)
#    print "vars: ", variables
#    print "prods: ", prods

    max_vars = reduction_possible(variables)
    new_prods = []
    no_mult = []
#    print "max_vars: ", max_vars
    if max_vars:
        for p in prods:
#            li = p.split(mult)
            li = get_products(p, format, mult)
            try:
                # See if we can find all variables in current product
                indices = [li.index(i) for i in max_vars]
            except:
                no_mult.append(p)
                continue
            for v in max_vars:
                li.remove(v)
            # If we've removed the last variable add 1.0
            if not li:
                li.append(str(1))
            p = mult.join(li)
            new_prods.append(p)
    else:
        # No reduction possible
        return expression
#    print "new_prods: ", new_prods
    new_prods = group_vars(add.join(new_prods), format)
    len_new_prods = len(get_products(new_prods, format, add))
#    print "new_prods: ", new_prods

#    l = mult.join([group % new_sum, max_var])
#    print l
#    no_mult += [l]
#    print no_mult
#    return "hej"
    # Recursively reduce sums with and without reduced variable

    if new_prods:
        # only pick the new string
        new_prods = reduce_operations(new_prods, format)
    if no_mult:
        # only pick the new string
        no_mult = [reduce_operations(add.join(no_mult), format)]

    g = new_prods
    if len_new_prods > 1:
        g = format["grouping"](new_prods)

    new_expression = add.join(no_mult + [mult.join([g, mult.join(max_vars)])])

    return new_expression

def expand_operations(expression, format):
    add   = format["add"](["", ""])
    mult  = format["multiply"](["", ""])
    div  = format["division"]
    group = format["grouping"]("")
    access = format["array access"]("")
    l = group[0]
    r = group[1]
#    print l
#    print r
    count = 0
#    print "\nexpr: ", expression
    # Check that we have the same number of left/right parenthesis in expression
    if not expression.count(l) == expression.count(r):
        raise RuntimeError, "Number of left/right parenthesis do not match"

    if expression.count(l) == 0:
        return group_vars(expression, format)

    prods = get_products(expression, format, add)
#    print "prods: ", prods
    new_prods = []
    if len(prods) > 1:
        for p in prods:
            new_prods.append(expand_operations(p, format))
        return group_vars(add.join(new_prods), format)

    if not len(prods) == 1:
        raise RuntimeError, "too many products"

    split_prods = get_products(prods[0], format, mult)
#    if not len(split_prods) > 1:
#        return group_vars(prods[0], format)

#    print "split prods: ", split_prods
    if len(split_prods) > 1:
        new_split = []
        for p in split_prods:
            if not l in p:
                new_split.append(p)
                split_prods.remove(p)
        if not split_prods:
            return group_vars(mult.join(new_split), format)
        else:
            fac = group_vars(mult.join(new_split), format)
#            print "fac: ", fac
            left = get_products(expand_operations(split_prods[0], format), format, add)
            right = get_products(expand_operations(mult.join(split_prods[1:]), format), format, add)
#            print "left: ", left
#            print "right: ", right
            new_mults = []
            for lp in left:
                for rp in right:
                    ent = [fac, lp, rp]
                    ent = [e for e in ent if e]
                    if ent:
                        new_mults.append(mult.join(ent))
#            print "new_mults: ", new_mults
            return group_vars(add.join(new_mults), format)

    if expression.count(l) > 1:
        nested = 0
        for c in expression:
            if c == l:
                nested += 1
            if c == r:
                nested -= 1
            if nested > 1:
                break
        if nested:
#            print "nested"
            new_expr = ""
            posr = 0
            count = 0
            par = False
            for c in expression:
                if c == l:
                    par = True
                    count += 1
                if c == r:
                    count -= 1
                if count == 0 and par:
                    break
                posr += 1
            posl = expression.index(l)
            inner = expression[posl+1:posr]
#            print "inner: ", inner
            new_inner = expand_operations(inner, format)
#            if len(get_products(new_inner, format, add)) == 1:
#                new_inner = new_inner[1:-1]
#            print "new_inner: ", new_inner
            # If nothing happened return same
            if inner == new_inner:
                return expression
            else:
                new_expr = expand_operations(expression.replace(inner, new_inner, 1), format)
                return group_vars(new_expr, format)

#    return expression
    if split_prods[0][0] == l and split_prods[0][-1] == r: 
        inner = split_prods[0][1:-1]
    else:
        inner = split_prods[0]
#    print "INNER: ", inner
    return inner
#    new_prods = []
#    ml = mult.join(["", l])
#    rm = mult.join([r,""])
#    for p in prods:
#        if r and l in p:
#            fac = []
#            if ml in p:
#                l_fac = get_products(p.split(ml)[0], format, add)
#                if len(l_fac) > 1:
#                    raise RuntimeError, "too many products"
#                fac += l_fac
##                print "l_fac: ", l_fac
#            if rm in p:
#                r_fac = get_products(p.split(rm)[-1], format, add)
#                if len(r_fac) > 1:
#                    raise RuntimeError, "too many products"
#                fac += r_fac
##                print "r_fac: ", r_fac
#            inner = get_products(p.split(r)[0].split(l)[-1], format, add)
##            print inner
#            new_p = add.join(inner)
#            print "fac: ", fac
#            if fac:
#                new_p = add.join([mult.join(fac + [i]) for i in inner])
#            new_prods.append(new_p)
#            print "new_p: ", new_p
#        else:
#            new_prods.append(p)
#    print "new_prods: ", new_prods
##    print "new_prods: ", add.join(new_prods)
##    print "new_prods: ", group_vars(add.join(new_prods), format)
#    new_expr = group_vars(add.join(new_prods), format)

##    new_expr = add.join(new_prods)

#    return new_expr

if __name__ == "__main__":

    simple_format = {
        "add": lambda v: " + ".join(v),
        "subtract": lambda v: " - ".join(v),
        "multiply": lambda v: "*".join(v),
        "division": "/",
        "grouping": lambda v: "(%s)" % v,
        "absolute value": lambda v: "std::abs(%s)" % v,
        "sqrt": lambda v: "std::sqrt(%s)" % v,
        "array access": lambda i: "[%s]" %(i)}

    expr0 = "(x*2) + y*x + x*x"     # should become:   (2 + y + x)*x     (5 --> 3)
    expr1 = "x + y*x + x*x"       # should become:   x + (y + x)*x     (4 --> 3)
    expr2 = "2*x*y + y*x*z + x*x" # should become:   ((2 + z)*y + x)*x (7 --> 4)
    expr3 = "2*x*y + y*x*3 + x*x" # should become:   (5*y + x)*x (7 --> 3)
    expr4 = "x*y + y*x + x*x"     # should become:   (2*y + x)*x (7 --> 3)
    expr5 = "y*x + x*y + x*y*z"  # should become:   (2 + z)*x*y (6 --> 3)
    expr6 = "P_t2_p0_s0_s0[ip][ind00[q]]*P_t2_p0_s0_s0[ip][ind00[w]]*W0[ip]*G_0 + P_t2_p0_s0_s0[ip][ind00[q]]*P_t2_p0_s0_s0[ip][ind00[w]]*W1[ip]*G_0000 + P_t2_p0_s0_s0[ip][ind00[q]]*P_t2_p0_s0_s0[ip][ind00[w]]*W2[ip]*G_0000"
    expr7 = "(2 + (7.0*(x + z)*y*y + 3.0)*z*y + (5) + x)*x" # should become:   2*x + x*x*y*z + 3*z*x + x*x (6 --> 10)
    expr8 = "P_t8_p1_s0[ip][nzc1[j]]*P_t11_p0_s1[ip][nzc0[k]]*W0[ip]*w[0][4 + ip]*Jinv_00*Jinv_11*det + P_t8_p1_s0[ip][nzc1[j]]*P_t11_p0_s1[ip][nzc0[k]]*W0[ip]*w[0][8 + ip]*Jinv_00*Jinv_10*det + P_t8_p1_s0[ip][nzc1[j]]*P_t11_p0_s1[ip][nzc0[k]]*W0[ip]*w[0][28 + ip]*Jinv_01*Jinv_11*det + P_t8_p1_s0[ip][nzc1[j]]*P_t11_p0_s1[ip][nzc0[k]]*W0[ip]*w[0][32 + ip]*Jinv_01*Jinv_10*det"
    expr9 = "P_t2_p1[ip][nzc1[j]]*P_t2_s0_a0[ip][nzc2[0]]*P_t2_p1[ip][nzc1[k]]*P_t2_s0_a0[ip][nzc2[0]]*W0[ip]*w[0][nzc2[0] + 3]*w[0][nzc2[0]]*det + P_t2_p1[ip][nzc1[j]]*P_t0_s1_a1[ip][nzc0[0]]*P_t2_p1[ip][nzc1[k]]*P_t0_s1_a1[ip][nzc0[0]]*W0[ip]*w[0][nzc0[0] + 3]*w[0][nzc0[0]]*det"

    expr10 = "P_t0_p1_s0_s0[ip][j]*P_t0_p1_s1_s0[ip][nzc0[k]]*W0[ip]*(J_00*Jinv_00*J_01*Jinv_00 + J_10*Jinv_00*J_11*Jinv_00 + J_00*Jinv_01*J_01*Jinv_01 + J_10*Jinv_01*J_11*Jinv_01)*det*1.0/(detJ*detJ)"# + P_t0_p1_s1_s1[ip][j]*P_t0_p1_s1_s0[ip][nzc0[k]]*W0[ip]*(J_01*Jinv_10*J_01*Jinv_00 + J_11*Jinv_10*J_11*Jinv_00 + J_01*Jinv_11*J_01*Jinv_01 + #J_11*Jinv_11*J_11*Jinv_01)*det*1.0/(detJ*detJ)"
    expr11 = "3*x + 3*(x + y)*(x + y)*(3)" # 3x + 9xx + 18xy + 9yy
#    expr12 = "3 + 2*(y + z)"# + 5*(x + y)*3 + 4" # 6x + 3y
#    expr12 = "y + z"# + 5*(x + y)*3 + 4" # 6x + 3y
#    expr12 = "2*((y + z) + 5) + 4" # 6x + 3y
#    expr12 = "2*x*y + 4*y*x + 5" # 6x + 3y
#    expr12 = "((x + y))*1.0/(4*y*x + 5)" # 6x + 3y
#    expr12 = "q*(2*x + (x + y)*(4*y + 5))" # 6x + 3y
#    expr12 = "(z + y)*1.0/(2*y + 5)" # 6x + 3y
#    expr12 = "w[i + 1]*1.0/(2*z + y) + 1.0/(2*y + 5)*w[2*i + 6]" # 6x + 3y
#    expr12 = "z*x + y + y" # 6x + 3y
#    expr12 = "(((y + z) + 5)) + 4" # 6x + 3y
#    expr12 = "P_t0_s0[ip][r]*P_t0_s0[ip][s]*P_t0_p0_s0[ip][nzc1[j]]*P_t0_p0_s0[ip][nzc1[k]]*W0[ip]*std::sqrt((1.0/std::abs((1.0/w[0][r]))))*std::sqrt(w[1][s])*(Jinv_00*Jinv_00 + Jinv_01*Jinv_01)*det"

    expr12 = "P_t0_p1_a0_s2[ip][nzc8[j]]*P_t2_p1_s1_s1[ip][nzc2[k]]*W1[ip]*Jinv_21*Jinv_10*det + P_t0_p1_a0_s2[ip][nzc8[j]]*P_t2_p1_s1_s1[ip][nzc2[k]]*W1[ip]*Jinv_21*Jinv_10*det"
#    expr12 = "((1.0/w[0][r]))"
#    print
#    print "reduce_operations(%s) --> (2 + y + x)*x" %expr0
#    res = reduce_operations(expr0, simple_format)
#    print "output:\n%d operations --> %d operations\nnew expression: %s" %\
#    (operation_count(expr0, simple_format), operation_count(res, simple_format), res)
#    print
#    print "reduce_operations(%s) --> x + (y + x)*x" %expr1
#    res = reduce_operations(expr1, simple_format)
#    print "output:\n%d operations --> %d operations\nnew expression: %s" %\
#    (operation_count(expr1, simple_format), operation_count(res, simple_format), res)
#    print
#    print "reduce_operations(%s) --> ((2 + z)*y + x)*x" %expr2
#    res = reduce_operations(expr2, simple_format)
#    print "output:\n%d operations --> %d operations\nnew expression: %s" %\
#    (operation_count(expr2, simple_format), operation_count(res, simple_format), res)
#    print
#    print "reduce_operations(%s) --> (5*y + x)*x" %expr3
#    res = reduce_operations(expr3, simple_format)
#    print "output:\n%d operations --> %d operations\nnew expression: %s" %\
#    (operation_count(expr3, simple_format), operation_count(res, simple_format), res)
#    print
#    print "reduce_operations(%s) --> (2*y + x)*x" %expr4
#    res = reduce_operations(expr4, simple_format)
#    print "output:\n%d operations --> %d operations\nnew expression: %s" %\
#    (operation_count(expr4, simple_format), operation_count(res, simple_format), res)
#    print
#    print "reduce_operations(%s) --> (2 + z)*y*x" %expr5
#    res = reduce_operations(expr5, simple_format)
#    print "output:\n%d operations --> %d operations\nnew expression: %s" %\
#    (operation_count(expr5, simple_format), operation_count(res, simple_format), res)
#    print
#    print "reduce_operations(%s) --> manual verification" %expr6
#    res = reduce_operations(expr6, simple_format)
#    print "output:\n%d operations --> %d operations\nnew expression: %s" %\
#    (operation_count(expr6, simple_format), operation_count(res, simple_format), res)
#    print
#    print expr7
#    print operation_count(expr7, simple_format)
#    print 
#    res_e = expand_operations(expr7, simple_format)
#    print res_e
#    print operation_count(res_e, simple_format)
#    print get_products("(2 + y)*(2 + 5 + 9 + 7)", simple_format, "*")
#    print get_products("x*y*1.0/(z + 5)", simple_format, "()")
#
#    print get_all_variables("x*y*1.0/(z + 5)", simple_format)

#    res_r = reduce_operations(res_e, simple_format)
#    print res_r
#    print operation_count(res_r, simple_format)
#    print 
#    print "reduce_operations(%s) --> manual verification" %expr8
#    res = reduce_operations(expr8, simple_format)
#    print "output:\n%d operations --> %d operations\nnew expression: %s" %\
#    (operation_count(expr8, simple_format), operation_count(res, simple_format), res)
#    print

    print expr12
    print
    expr12 = expand_operations(expr12, simple_format)
    print "expr12: ", expr12
    print
#    print group_vars("x + 1*x", simple_format)
    print reduce_operations(expr12, simple_format)
#    geo = {}
#    print get_geo_terms(expr8, geo, simple_format)
#    print geo

#    print collect_floats("2*4.0*x + 5*x", simple_format)
#    print reduce_operations("2*4.0*x + 5*x", simple_format)
#    res_r = reduce_operations(expr7, simple_format)
#    print res_r
#    print "expand_operations(%s) --> %s" %(res_r, expr0)
#    print "ops reduced: %s is %d" %(res_r, operation_count(res_r, simple_format))
#    print "ops expanded: %s is %d" %(res_e, operation_count(res_e, simple_format))
#    res = reduce_operations(expr0, simple_format)
#    print "output:\n%d operations --> %d operations\nnew expression: %s" %\
#    (operation_count(expr0, simple_format), operation_count(res, simple_format), res)
    print











