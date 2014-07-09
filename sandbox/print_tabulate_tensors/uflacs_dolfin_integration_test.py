from dolfin import *

output = {}
for representation in ("uflacs", "quadrature"):
    parameters["form_compiler"]["representation"] = representation
    print "Using form compiler representation:", representation
    output[representation] = {}
    for n in (1, 2, 4):
        print "n =", n
        out = []
        output[representation][n] = out

        mesh = UnitSquareMesh(n,3)
        V = FunctionSpace(mesh, "CG", 1)
        c = Constant(2.3)
        f = Function(V)
        e = Expression("x[0]+x[1]")
        f.interpolate(e)
        u = TrialFunction(V)
        v = TestFunction(V)

        out.append( assemble(c*dx, mesh=mesh) )
        out.append( assemble(f*dx) )
        out.append( assemble(f**2*dx) )
        out.append( assemble(c*f*dx) )
        out.append( assemble(grad(f)**2*dx) ) # Fails!

        out.append( assemble(v*dx).norm('l2') )
        out.append( assemble(c*v*dx).norm('l2') )
        out.append( assemble(v.dx(0)*dx).norm('l2') )

        out.append( assemble(u*v*dx).norm('frobenius') )
        out.append( assemble(dot(grad(u),grad(v))*dx).norm('frobenius') )

	m = len(out)


for i in xrange(m):
    print
    print "i =", i
    for n in (1, 2, 4):
        u = output["uflacs"][n][i]
        q = output["quadrature"][n][i]
        d = q-u
        a = (abs(q)+abs(u))/2
        r = 2*abs(d) / a
        if a > 0.0 and r > 1e-8:
            print n, ("%10.2e"%(d)), ("%10.2e"%(r)), ("%10.2e"%u), ("%10.2e"%q)
        else:
            print n, "ok!"

