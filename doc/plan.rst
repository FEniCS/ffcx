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

- Generate fully usable dolfin Expressions and test with PyDOLFIN

- Make jit test framework to simplify generate-compile-build-run-test cycle

- Start generating dune code

