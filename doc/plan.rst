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

uflacs compile_pdelab <files>
-----------------------------
TODO: Compile forms into PDELAB/DUNE C++ code.


Backlog
=======

Make a usable compiler
----------------------

- Generate fully usable dolfin Expressions and test with PyDOLFIN

- Generate forms from sfc with uflacs generic compiler and test with PyDOLFIN

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

Testing
-------

- Make a faster test framework for compiled code not using swig,
  to simplify and speed up generate-compile-build-run-test cycle

Plugins
-------

Potential plugins, probably won't implement all these.

- Implement an itk compiler plugin

- Implement a vtk compiler plugin

- Implement a vmtk compiler plugin

- Implement a dune compiler plugin

- Implement a diffpack compiler plugin

- Implement a deal.II compiler plugin
