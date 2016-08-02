.. title:: FFC


=============================
FFC: The FEniCS Form Compiler
=============================

FFC is a compiler for finite element variational forms. From a
high-level description of the form, it generates efficient low-level
C++ code that can be used to assemble the corresponding discrete
operator (tensor). In particular, a bilinear form may be assembled
into a matrix and a linear form may be assembled into a vector.  FFC
may be used either from the command line (by invoking the ``ffc``
command) or as a Python module (``import ffc``).

FFC is part of the FEniCS Project.

For more information, visit http://www.fenicsproject.org


Documentation
=============

.. toctree::
   :titlesonly:
   :maxdepth: 1

   installation
   manual
   API reference (FFC) <api-doc/ffc>
   API reference (UFLACS) <api-doc/uflacs>
   releases

[FIXME: These links don't belong here, should go under API reference somehow.]

* :ref:`genindex`
* :ref:`modindex`
