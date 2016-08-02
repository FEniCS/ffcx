.. title:: Installation


============
Installation
============

FFC is normally installed as part of an installation of FEniCS.
If you are using FFC as part of the FEniCS software suite, it
is recommended that you follow the
`installation instructions for FEniCS
<https://fenics.readthedocs.io/en/latest/>`__.

To install FFC itself, read on below for a list of requirements
and installation instructions.

Requirements and dependencies
=============================

FFC requires Python version 2.7 or later and depends on the
following Python packages:

* NumPy
* SymPy
* six

These packages will be automatically installed as part of the
installation of FFC, if not already present on your system.

Installation instructions
=========================

To install FFC, download the source code from the
`FFC Bitbucket repository
<https://bitbucket.org/fenics-project/ffc>`__,
and run the following command:

.. code-block:: console

    pip install .

To install to a specific location, add the ``--prefix`` flag
to the installation command:

.. code-block:: console

    pip install --prefix=<some directory> .
