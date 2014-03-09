#!/usr/bin/env python

import sys
sys.path.insert(0, "..")
import uflacs
print uflacs.__file__

from uflacs import *


from ufl import *
from ufl.classes import GeometricQuantity
from ufl.algorithms import replace, change_to_reference_grad

if 0:
    domain0 = Domain(triangle)
    V0 = VectorElement("CG", domain0, 1)
    x = Coefficient(V0)
    domain = Domain(x)
    # This currently seems to fail in ufl, thought it was supposed to work now...
else:
    domain = Domain(triangle)

V = FiniteElement("CG", domain, 2)
dx = Measure("cell", domain)

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
if case == 4:
    M = grad(f)[0]*dx
    L = grad(v)[0]*dx
    a = dot(2*grad(u),f*grad(v))*dx

forms = [M, L, a]

from uflacs.representation.compute_expr_ir import *

for form in forms:
    print '/'*80

    expr = form.integrals()[0].integrand()
    #print "Initial expression"
    #print str(expr)

    fd = form.compute_form_data()
    expr = fd.preprocessed_form.integrals()[0].integrand()
    #print "First apply ufl preprocessing to form"
    #print str(expr)


    # FIXME: Apply these transformations at a suitable place in uflacs compiler
    expr = replace(expr, fd.function_replace_map)
    #print "Then function replace map"
    #print str(expr)

    expr = change_to_reference_grad(expr)
    print "And change to local grad"
    print str(expr)


    #print "Build list based graph representation of scalar subexpressions"
    expressions = [expr]

    e2i, V, target_variables, modified_terminals = build_scalar_graph(expressions)

    if 0:
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
    argument_factorization, modified_arguments, V, target_variables, dependencies = \
        compute_argument_factorization(V, target_variables, dependencies)
    # modified_arguments = [v, ...] where each v is a modified argument
    # V = [v, ...] where each v is argument independent
    # factorization = { (i,...): j } where (i,...) are indices into modified_arguments and j is an index into factorized_vertices
    if 0:
        print
        print '\nargument_factorization'
        print format_mapping(argument_factorization)
        print '\nargument factors'
        print format_enumerated_sequence(modified_arguments)
        print '\nV'
        print format_enumerated_sequence(V)
        print '\ntarget_variables'
        print target_variables
        print '\ndependencies'
        print dependencies
        print

    # Count the number of dependencies every subexpr has
    depcount = compute_dependency_count(dependencies)

    # Build the 'inverse' of the sparse dependency matrix
    inverse_dependencies = invert_dependencies(dependencies, depcount)

    print "Mark subexpressions of V that are actually needed for final result"
    active, num_active = mark_active(dependencies, target_variables)

    print "Build set of modified_terminal indices into factorized_vertices"
    modified_terminal_indices = [i for i,v in enumerate(V)
                                 if is_modified_terminal(v)]

    print "Build piecewise/varying markers for factorized_vertices"
    spatially_dependent_terminal_indices = [i for i in modified_terminal_indices
                                   if not V[i].is_cellwise_constant()]
    spatially_dependent_indices, num_spatial = mark_image(inverse_dependencies, spatially_dependent_terminal_indices)

    print
    print "Vertices:"
    print format_enumerated_sequence(V)
    print
    print "Active:", num_active, len(V)
    print format_enumerated_sequence(active)
    print
    print "Modified terminals:"
    print modified_terminal_indices
    print
    print "Spatially dependent:"
    print spatially_dependent_terminal_indices
    print
    print spatially_dependent_indices
    print


    # This is returned from ir compiler:
    expr_ir = {}
    expr_ir["argument_factorization"] = argument_factorization
    expr_ir["modified_arguments"] = modified_arguments
    expr_ir["V"] = V
    expr_ir["target_variables"] = target_variables
    expr_ir["active"] = active
    expr_ir["dependencies"] = dependencies
    expr_ir["inverse_dependencies"] = inverse_dependencies
    expr_ir["modified_terminal_indices"] = modified_terminal_indices
    expr_ir["spatially_dependent_indices"] = spatially_dependent_indices



    # ... Tables enter here
    from uflacs.analysis.modified_terminals import analyse_modified_terminal2
    #print "Build modified_argument_tables = { modified_arguments_index: (uname,b,e) }"
    modified_argument_tables = {}
    for i, a in enumerate(expr_ir["modified_arguments"]):
        mt = analyse_modified_terminal2(a)
        (uname, b, e) = (str(a), mt.terminal.number(), mt.terminal.number()+3) # FIXME
        modified_argument_tables[i] = (uname, b, e)

    #print "Build modified_terminal_tables = { factorized_vertices_index: (uname,b,e) } for coeffs,jacobian"
    #gtypes = (Jacobian, JacobianInverse, JacobianDeterminant, SpatialCoordinates, LocalCoordinates)
    modified_terminal_tables = {}
    for i in expr_ir["modified_terminal_indices"]:
        v = expr_ir["V"][i]
        mt = analyse_modified_terminal2(v)
        element = None
        if isinstance(mt.terminal, Coefficient):
            element = mt.terminal.element()
        elif isinstance(mt.terminal, GeometricQuantity): # TODO: Refine this check with gtypes?
            D = mt.terminal.domain()
            coordinates = D.coordinates()
            if coordinates:
                element = coordinates.element()
        if element is not None:
            (uname, b, e) = (str(v), 0, 0) # FIXME
            modified_terminal_tables[i] = (uname, b, e)

    print
    print "modified_argument_tables:"
    print format_mapping(modified_argument_tables)
    print "modified_terminal_tables:"
    print format_mapping(modified_terminal_tables)

    # ... Code generation starts here

    #print "TODO: Generate code for defining tables referenced by"\
    #      " modified_argument_tables and modified_terminal_tables"

    #print "TODO: Generate code for loop nests in tabulate_tensor with"\
    #      " blocks of A[(i0+b0)*n1+(i1+b1)] += f*v0[i0]*v1[i1]"

    # Find the blocks to build:
    arguments_and_factors = sorted(expr_ir["argument_factorization"].items(), key=lambda x: x[0])
    for args, factor_index in arguments_and_factors:
        dofblock = ()
        argunames = ()
        for arg in args:
            uname, b, e = modified_argument_tables[arg]
            dofblock = dofblock + ((b,e),)
            argunames = argunames + (uname,)
        print "\nTODO: Emit code for {dofblock}, {argunames}, {factor_index}".format(**locals())


    #print "TODO: Generate code for cellwise and spatial partitions of factorized_vertices"

    #print "TODO: Generate code for coefficients using modified_terminal_tables"

    #print "TODO: Generate code for geometry"

    """How to generate each modified terminal type:

    - Argument: table name and dofrange says it all, generate directly from table name and index names

    - Coefficient[degree==0]: use dofname directly
    - Coefficient[degree>0]: generate (inline) dot product between dofnames and table

    - Jacobian: generate (inline) dot product between "vertex dofnames" and ReferenceGrad

    - JacobianDeterminant: generate code snippet for now
    - JacobianInverse: generate code snippet for now

    - SpatialCoordinate: generate (inline) dot product between "vertex dofnames" and table
    - LocalCoordinate: Need FacetLocalCoordinate as well. Make all Facet* subclass FacetGeometricQuantity?
    """

    print '\\'*80

"""How to approximate tensor representation as a special case:

 - Choose to compute factors first for all quadrature points:
   for(iq)
   {
     ...;
     F[iq][k] = ...;
   }

 - Then implement the tensor contractions smartly:
   A[i0*n1 + i1] = sum(
       V0[k][iq][i0-b0[k]] * V1[k][iq][i1-b1[k]] * F[k][iq] * weights[iq]
       for i0 in [b0[k], e0[k])
       for i1 in [b1[k], e1[k])
       for iq in [0,nq)
       for k in [0,nterms)
       ) * detJ;

 - We can always precompute product of arguments and weights at the cost of large tables:
   - Precompute: T[i0][i1][k][iq] = V0[k][iq][i0-b0[k]] * V1[k][iq][i1-b1[k]] * weights[iq]
   - Compute: F[k][iq] = ...; with arbitrary code (not limited like tensor representation)
   - Compute: A[i0][j0] = sum(T[i0][j0][k][iq]*F[k][iq] for k in terms for iq in [0,nq))

 - If F[k][iq] = F[k], we can precompute the sum of T over iq:
   - T[k][i0][i1] = sum(V0[k][iq][i0-b0[k]] * V1[k][iq][i1-b1[k]] * weights[iq] for iq in [0,nq))
   - Compute: F[k] = ...; with arbitrary code (not limited like tensor representation)
   - Compute: A[i0][j0] = sum(T[i0][j0][k]*F[k] for k in terms)

 - If F[k][iq] = F[k] = F, we can do fully precomputable mass matrix:
   - T[i0][i1] = sum(V0[k][iq][i0-b0[k]] * V1[k][iq][i1-b1[k]] * weights[iq] for iq in [0,nq) for k in terms)
   - Compute: F = ...; with arbitrary code (not limited like tensor representation)
   - Compute: A[i0][j0] = T[i0][j0] * F
"""
