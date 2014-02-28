#!/usr/bin/env python

import sys
sys.path.insert(0, "..")
import uflacs
print uflacs.__file__

from uflacs import *


from ufl import *
from ufl.algorithms import replace, change_to_local_grad

domain = Domain(triangle)
dx = Measure("cell", domain)
V = FiniteElement("CG", domain, 1)

u = TrialFunction(V)
v = TestFunction(V)
c = Constant(domain)
f = Coefficient(V)

case = int(sys.argv[1])

if case == 1:
    M = 1*dx
    L = v*dx
    a = u*v*dx
if case == 2:
    M = f*dx
    L = f*v*dx
    a = f*u*v*dx
if case == 3:
    M = grad(f)[0]*dx
    L = grad(v)[0]*dx
    a = dot(grad(u),grad(v))*dx

forms = [M, L, a]

from uflacs.generation.compiler import *

for form in forms:
    print '/'*80

    print "Initial expression"
    expr = form.integrals()[0].integrand()
    print str(expr)

    print "First apply ufl preprocessing to form"
    fd = form.compute_form_data()
    expr = fd.preprocessed_form.integrals()[0].integrand()
    print str(expr)

    print "Then function replace map"
    expr = replace(expr, fd.function_replace_map)
    print str(expr)

    print "And change to local grad"
    expr = change_to_local_grad(expr)
    print str(expr)

    #print "TODO: Build list based graph representation of scalar subexpressions"
    expressions = [expr]

    e2i, V, target_variables, modified_terminals = build_scalar_graph(expressions)

    print
    print "\nV:"
    print format_enumerated_sequence(V)
    print "\ne2i:"
    print format_mapping(e2i)
    print "\ntarget_variables:"
    print format_enumerated_sequence(target_variables)
    print "\nmodified_terminals:"
    print format_enumerated_sequence(modified_terminals)
    print

    dependencies = compute_dependencies(e2i, V)
    print '\ndependencies:'
    print format_enumerated_sequence(dependencies)

    print "Build factorization"
    # AV, FV, IM
    argument_factors, factorized_vertices, argument_factorization, target_variables, dependencies = \
        compute_argument_factorization(V, target_variables, dependencies)
    # argument_factors = [v, ...] where each v is a modified argument
    # factorized_vertices = [v, ...] where each v is argument independent
    # factorization = { (i,...): j } where (i,...) are indices into argument_factors and j is an index into factorized_vertices
    print
    print '\nargument factors'
    print format_enumerated_sequence(argument_factors)
    print '\nfactorized_vertices'
    print format_enumerated_sequence(factorized_vertices)
    print '\nargument_factorization'
    print format_mapping(argument_factorization)
    print '\ntarget_variables'
    print target_variables
    print '\ndependencies'
    print dependencies
    print


    # This is crap:
    # Rebuild some graphs from factorization
    #V, e2i, dependencies = rebuild_scalar_graph_from_factorization(
    #    argument_factors, factorized_vertices, argument_factorization)
    # TODO: target_variables for non-scalar or multiple expressions
    #target_variables = [len(V)-1]


    #print "TODO: Build piecewise/varying markers for factorized_vertices"

    #print "TODO: Build set of modified_terminal indices into factorized_vertices"

    # ... Tables enter here

    #print "TODO: Build modified_argument_tables = { argument_factors_index: (uname,b,e) }

    #print "TODO: Build modified_terminal_tables = { factorized_vertices_index: (uname,b,e) } for args,coeffs,jacobian"

    #print "TODO: Generate code for defining tables referenced by modified_argument_tables and modified_terminal_tables"

    #print "TODO: Generate code for loop nests in tabulate_tensor with blocks of A[(i0+b0)*n1+(i1+b1)] += f*v0[i0]*v1[i1]"

    #print "TODO: Generate code for cellwise and spatial partitions of factorized_vertices"

    #print "TODO: Generate code for coefficients using modified_terminal_tables"

    #print "TODO: Generate code for geometry"

    print '\\'*80
