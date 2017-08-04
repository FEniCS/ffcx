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
* six

FFC also depends on the following FEniCS Python packages:

* FIAT
* UFL
* dijitso

These packages will be automatically installed as part of the
installation of FFC, if not already present on your system.

.. _tsfc_requirements:

TSFC requirements
-----------------

To use experimental ``tsfc`` representation, additional
dependencies are needed:

* `TSFC <https://github.com/blechta/tsfc>`_ [1]_
* `COFFEE <https://github.com/blechta/COFFEE>`_ [1]_
* `FInAT <https://github.com/blechta/FInAT>`_ [1]_

and in turn their additional dependencies:

* singledispatch [2]_
* networkx [2]_
* PuLP [2]_

.. note:: TSFC requirements are not installed in FEniCS Docker
    images by default yet but they can be easilly installed
    on demand::

        docker pull quay.io/fenicsproject/dev:latest
        docker run -ti --rm quay.io/fenicsproject/dev:latest
        sudo apt-get update && sudo apt-get -y install glpk-utils && \
          pip2 install --prefix=${FENICS_PREFIX} --no-cache-dir \
          git+https://github.com/blechta/tsfc.git \
          git+https://github.com/blechta/COFFEE.git \
          git+https://github.com/blechta/FInAT.git \
          singledispatch networkx pulp && \
          pip3 install --prefix=${FENICS_PREFIX} --no-cache-dir \
          git+https://github.com/blechta/tsfc.git \
          git+https://github.com/blechta/COFFEE.git \
          git+https://github.com/blechta/FInAT.git \
          singledispatch networkx pulp && \
          sudo apt-get clean && \
          sudo rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

    The first two commands (or their modification, or
    ``fenicsproject`` helper script) are to be run on a host,
    while the last command, to be run in the container, actually
    installs all the TSFC requirements. For further reading,
    see `FEniCS Docker reference
    <https://fenics-containers.readthedocs.io/>`_.

.. [1] These are forks of the original packages tested to be
   compatible with FFC and updated frequently from upstream.

.. [2] Pip-installable.

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
