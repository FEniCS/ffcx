FFC/UFC element factory module (ffc_test_factory)
=================================================

This module provides a factory for FFC generated code objects. UFC
code is generated, and then wrapped using pybind11.


Creating the UFC and factory interface code
-------------------------------------------

Running the Python script

   python3 ./generate_factory.py

creates a pybind11 source file for building the factory module.


Building and installing the factory module
------------------------------------------

Requirement: `ufc_wrappers` module and `pybind11`.

To build and install the factory:

    pip3 install .


Extending the factory
---------------------

The Python script `generate_factory.py` can be extended to add more
UFC objects to the factory.
