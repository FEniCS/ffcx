UFLAC - UFL Analyser and Compiler
=================================


Commands
========

uflacs load <files>
-------------------
Attempt loading files.
Report any errors.
Report information about the file contents.

uflacs repr <files>
-------------------
Load files and print repr of all loaded objects.

uflacs str <files>
------------------
Load files and print str of all loaded objects.

uflacs tree <files>
-------------------
Load files and print tree formatting of all loaded objects.

uflacs analyse <files>
----------------------
Attempt preprocessing of files.
Report any errors.
Report information about the file contents and result of preprocessing.

uflacs latex <files>
--------------------
Compile files into a latex document.
TODO: Render elements.
TODO: Render expressions.
TODO: Improve form rendering.

uflacs graphviz <files>
-----------------------
TODO: Compile files into graphviz .dot graph files for visualization.

uflacs compile <files>
----------------------
TODO: Compile forms into half finished C++ code.

uflacs compile_dolfin <files>
-----------------------------
TODO: Compile expressions into DOLFIN C++ code.


Backlog
=======

Testing
-------

- Use cmake to get ufc paths in test framework

- Use cmake to get dolfin paths in test framework and link with dolfin for testing Expression etc?

- Make test framework for compiled code more robust w.r.t. code changes and rerunning

Make a usable ufl-to-C++ compiler
---------------------------------

- Generate code following ufc conventions for all geometric ufl quantities

- Generate dolfin Expressions and test with dolfin

- Generate functionals from ffc and test with dolfin

- Generate linear forms from ffc

- Generate bilinear forms from ffc

- Generate n-linear forms from ffc


Improve compiler algorithms
---------------------------

- Rewrite integrand terms on the form "((factor) * D(u,n)) * D(u,m)"
  to optimize loops

- Rewrite integrand terms on the form "((factor) * D(u,n)) * D(u,m)"
  to optimize loops

- Improve compiler plugin interface. Make default version configurable
  so only minimal extra code is necessary.

Improve cache score heuristics
..............................

- Let literal*f count less than f*g, e.g. by using cost of operands in
  heuristic, but keep in mind that exp(f,literal) is still as costly
  as exp(f,g)!

- If expr2 in partition2 depends on expr1 in partition1, with
  partition1 < partition2, let that dependency count more.  Maybe use
  a score like (partition2 - partition1 + 1), or (partition2 -
  partition1 + 1)**2, or c1*(partition2 - partition1 + c2)**c3.

- Iterative algorithm? Pick x% highest scores, update rest of the
  scores, repeat until end criteria.

- Expose tuning parameters to external compiler interface

Improve register allocation algorithm
.....................................

- Store and allocate registers separately for each partition,
  i.e. declare "double s%(p)d[%(n)d];" % (p,num_registers[p])
  at the beginning of each partition

- Allow reuse of registers within partitions.

Plugins
-------

Potential plugins, probably won't implement all these.

- Implement a vtk compiler plugin

- Implement an itk compiler plugin

- Implement a vmtk compiler plugin

- Implement a diffpack compiler plugin

- Implement a dune/pdelab compiler plugin

- Implement a deal.II compiler plugin

