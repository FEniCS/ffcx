# Copyright (c) 2021-2022 Chris Richardson, Garth Wells, Matthew Scroggs
# FEniCS Project
# SPDX-License-Identifier: MIT

import pytest
import pkg_resources
import os


def test_version_numbering():
    py_version = pkg_resources.get_distribution("fenics-ffcx").version
    cpp_version = py_version.replace('dev', '')

    path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "..")
    if not os.path.isfile(os.path.join(path, "LICENSE")):
        pytest.skip("This test can only be run from the source directory.")

    print("Checking version numbers in cmake/CMakeLists.txt")
    found = False
    with open(os.path.join(path, "cmake/CMakeLists.txt")) as f:
        for line in f:
            if "ufcx version" in line.lower():
                found = True
                v = line.lower().split("ufcx version")[1].strip().split(" ")[0].strip()
                assert v == cpp_version
    assert found

    print("Checking version numbers in setup.cfg")
    found = False
    with open(os.path.join(path, "setup.cfg")) as f:
        for line in f:
            if "version = " in line:
                found = True
                v = line.split("version = ")[1].strip()
                assert v == py_version
    assert found

    print("Checking version numbers in ffcx/codegeneration/ufcx.h")
    with open(os.path.join(path, "ffcx/codegeneration/ufcx.h")) as f:
        content = f.read()
    ufcx_version = ".".join([
        content.split(f"#define UFCX_VERSION_{i} ")[1].split("\n")[0].strip()
        for i in ["MAJOR", "MINOR", "MAINTENANCE", "RELEASE"]
    ])
    assert ufcx_version == cpp_version
