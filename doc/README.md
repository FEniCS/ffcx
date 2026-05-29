# Building the FFCx Python documentation

To build the documentation:

1. Install FFCx using the ``docs`` optional dependency set, e.g.

       python -m pip install .[docs]

2. Run in this directory:
 
       python -m sphinx -W -b html source/ build/html/
