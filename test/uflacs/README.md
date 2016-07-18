# Test structure for UFLACS

- unit/

    unit tests of internal components of UFLACS

- crosslanguage/

    unit tests which produce C++ tests of generated code which is then
    executed by Google Test

- system/

    tests that use external software with uflacs, in particular
    integration with DOLFIN


## Running examples

    cd test/
    py.test
    py.test unit
    py.test system
    py.test crosslanguage
