#!/usr/bin/env python

def operation_count(expression, format):
    # Cheat to get character for add and multiply
    add   = format["add"](["", ""])
    mult  = format["multiply"](["", ""])
    adds = len(get_products(expression, format)) - 1
    return expression.count(mult) + adds

def get_products(expression, format):
    add   = format["add"](["", ""])
    mult  = format["multiply"](["", ""])
    access = format["array access"]("")
    l = access[0]
    r = access[1]
#    print expression
    new_prods = []
    incomplete = []
    prods = expression.split(add)

    while prods:
        p = prods[0]
#        print p
        pos_rl = p.rfind(l)
        pos_rr = p.rfind(r)
        left = False
        # Found both
        if not pos_rl == -1 and not pos_rr == -1:
            pos_l = p.find(l)
            pos_r = p.find(r)
            if pos_rr < pos_rl:
                left = True
                incomplete.append(p)
#                print "left side"
#                print p
            elif pos_r < pos_l:
                incomplete.append(p)
                if not left:
                    new_prods.append(add.join(incomplete))
#                    print "new_prods: ", new_prods
                    incomplete = []
#                    print "right side"
#                    print p
            else:
#                print "else: "
                new_prods.append(p)

        elif not pos_rl == -1 and pos_rr == -1:
            left = True
            incomplete.append(p)
#            print "left side"
#            print p
        elif pos_rl == -1 and not pos_rr == -1:
            incomplete.append(p)
            if not left:
                new_prods.append(add.join(incomplete))
#                print "new_prods: ", new_prods
                incomplete = []
#                print "right side"
#                print p

#            p = add.join([p, prods[i+1]])
#            new_prods.append(p)
#            # Increment i because we already took the next product
#            i += 1
        else:
#            print "else: "
            new_prods.append(p)

        prods.remove(p)

    return new_prods

def collect_floats(expression, format):
    add   = format["add"](["", ""])
    mult  = format["multiply"](["", ""])

    floats = 0.0
    new_prods = []
    for p in get_products(expression, format):
        try:
            floats += float(eval(p))
        except:
            new_prods.append(p)
    if floats:
        new_prods.append(str(floats))
    return add.join(new_prods)

def get_variables(expression, format):
    add   = format["add"](["", ""])
    mult  = format["multiply"](["", ""])

    prods = get_products(expression, format)
#    print "prods: ", prods
    variables = {}
    for i in range(len(prods)):
        p = prods[i]
        # only extract unique variables
        vrs = list(set(p.split(mult)))
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

    expression = expand_operations(expression, format)
    prods, variables = get_variables(expression, format)
#    print "vars: ", variables

    max_vars = reduction_possible(variables)
    new_prods = []
    no_mult = []
#    print "max_var: ", max_var
    if max_vars:
        for p in prods:
            li = p.split(mult)
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
#    print "new_sums: ", new_sums
    new_prods = collect_floats(add.join(new_prods), format)
    len_new_prods = len(get_products(new_prods, format))
#    new_sum = add.join(new_sums)
#    print no_mult
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
    group = format["grouping"]("")
    l = group[0]
    r = group[1]
#    print l
#    print r
    count = 0
    # Check that we have the same number of left/right parenthesis in expression
    if not expression.count(l) == expression.count(r):
        raise RuntimeError, "Number of left/right parenthesis do not match"

    pos_r = expression.find(r)
#    print "pos_r: ", pos_r

    if pos_r == -1:
        return expression

    # Get part of expression on the left side of ")" (includes the inner part)
    left0 = expression[0:pos_r]
    # Get factors from the RHS of ")"
    right = expression[pos_r+1 :]#.split(r)[0].split(add)[0]
    right_fac = expression[pos_r+1 :].split(r)[0].split(add)[0]
    pos_l = left0.rfind(l)
#    print "pos_l: ", pos_l
    # Get the inner part that needs to be factored
#    inner = expression[pos_l+1:pos_r]
    inner = left0.split(l)[-1]

    # Get the factors from the LHS of "("
    left = left0[0:pos_l]#.split(l)[-1].split(add)[-1]
    left_fac = left0[0:pos_l].split(l)[-1].split(add)[-1]
#    print "left: ", left
#    print "right: ", right
#    print "inner: ", inner
    factor = left_fac.split(mult) + right_fac.split(mult)
    factor = format["multiply"]([f for f in factor if not f == ""])
#    print "factor: ", factor
    # Multiply all terms in parenthesis with factor
    new_inner = ""
    if factor:
      new_inner = format["add"]([format["multiply"]([factor, s]) for s in inner.split(add)])
    else:
      new_inner = inner

    inner_replace = left_fac + format["grouping"](inner) + right_fac
#    print "new_inner", new_inner
#    print "inner_replace: ", inner_replace
#    print "expr find", expression.find(inner_replace)
    new_expr = expression.replace(inner_replace, new_inner)
    # Check if we still have some parenthesis left
    if new_expr.find(l) == -1:
        return new_expr
    else:
        return expand_operations(new_expr, format)
#    expression.replace(" ", "h")
#    print "new_expr: ", new_expr

#    return new_expr

if __name__ == "__main__":

    simple_format = {
        "add": lambda v: " + ".join(v),
        "subtract": lambda v: " - ".join(v),
        "multiply": lambda v: "*".join(v),
        "grouping": lambda v: "(%s)" % v,
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

    print
    print "reduce_operations(%s) --> (2 + y + x)*x" %expr0
    res = reduce_operations(expr0, simple_format)
    print "output:\n%d operations --> %d operations\nnew expression: %s" %\
    (operation_count(expr0, simple_format), operation_count(res, simple_format), res)
    print
    print "reduce_operations(%s) --> x + (y + x)*x" %expr1
    res = reduce_operations(expr1, simple_format)
    print "output:\n%d operations --> %d operations\nnew expression: %s" %\
    (operation_count(expr1, simple_format), operation_count(res, simple_format), res)
    print
    print "reduce_operations(%s) --> ((2 + z)*y + x)*x" %expr2
    res = reduce_operations(expr2, simple_format)
    print "output:\n%d operations --> %d operations\nnew expression: %s" %\
    (operation_count(expr2, simple_format), operation_count(res, simple_format), res)
    print
    print "reduce_operations(%s) --> (5*y + x)*x" %expr3
    res = reduce_operations(expr3, simple_format)
    print "output:\n%d operations --> %d operations\nnew expression: %s" %\
    (operation_count(expr3, simple_format), operation_count(res, simple_format), res)
    print
    print "reduce_operations(%s) --> (2*y + x)*x" %expr4
    res = reduce_operations(expr4, simple_format)
    print "output:\n%d operations --> %d operations\nnew expression: %s" %\
    (operation_count(expr4, simple_format), operation_count(res, simple_format), res)
    print
    print "reduce_operations(%s) --> (2 + z)*y*x" %expr5
    res = reduce_operations(expr5, simple_format)
    print "output:\n%d operations --> %d operations\nnew expression: %s" %\
    (operation_count(expr5, simple_format), operation_count(res, simple_format), res)
    print
    print "reduce_operations(%s) --> manual verification" %expr6
    res = reduce_operations(expr6, simple_format)
    print "output:\n%d operations --> %d operations\nnew expression: %s" %\
    (operation_count(expr6, simple_format), operation_count(res, simple_format), res)
    print
    print expr7
    print operation_count(expr7, simple_format)
    print 
    res_e = expand_operations(expr7, simple_format)
    print res_e
    print operation_count(res_e, simple_format)
    print 
    res_r = reduce_operations(res_e, simple_format)
    print res_r
    print operation_count(res_r, simple_format)
    print 
    print "reduce_operations(%s) --> manual verification" %expr8
    res = reduce_operations(expr8, simple_format)
    print "output:\n%d operations --> %d operations\nnew expression: %s" %\
    (operation_count(expr8, simple_format), operation_count(res, simple_format), res)
    print


#    print get_products(expr8, simple_format)

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











