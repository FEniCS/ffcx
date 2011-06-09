UFLAC - UFL Analyser and Compiler
=================================


Commands
========

uflacs load <files>
------------------
Attempt loading files.
Report any errors.
Report information about the file contents.

uflacs repr <files>
------------------
Load files and print repr of all loaded objects.

uflacs str <files>
-----------------
Load files and print str of all loaded objects.

uflacs latex <files>
-------------------
Compile files into a latex document.

uflacs dot <files>
-----------------
Compile files into graphviz .dot graph files for visualization.

uflacs analyse <files>
---------------------
Attempt preprocessing of files.
Report any errors.
Report information about the file contents and result of preprocessing.

uflacs compile <files>
---------------------
TODO Compile forms into C++ code templates.

uflacs compile-dolfin-expressions <files>
----------------------------------------
TODO Compile expressions into DOLFIN C++ code.

uflacs compile-dune <files>
--------------------------
TODO Compile forms into DUNE C++ code.


Backlog
=======

- Fetch code formatting utilities into repo

- Fill in some tests in test framework

- Make test framework to simplify generate-compile-build-run-test cycle

- Generate dolfin Expressions

- Start generating dune code

