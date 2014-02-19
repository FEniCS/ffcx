#!/bin/bash

RESULT=0

rm -f generated/*

rm -f python_tests.log
pushd ../test
#py.test 2>&1 | tee -a ../tests/python_tests.log
py.test &> ../tests/python_tests.log
if [ $? -ne 0 ]; then
    RESULT=1
    echo Python tests FAILED, see python_tests.log.
else
    echo Python tests PASSED.
fi
popd

rm -f cpp_tests.log
#make 2>&1 | tee -a cpp_tests.log
make &> cpp_tests.log
if [ $? -ne 0 ]; then
    RESULT=2
    echo C++ tests FAILED, see cpp_tests.log.
else
    echo C++ tests PASSED.
fi

exit $RESULT

