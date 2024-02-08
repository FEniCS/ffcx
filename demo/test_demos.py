"""Test demos."""

import os
import sys

import pytest

demo_dir = os.path.dirname(os.path.realpath(__file__))

ufl_files = []
for file in os.listdir(demo_dir):
    if file.endswith(".py") and not file == "test_demos.py":
        ufl_files.append(file[:-3])


@pytest.mark.parametrize("file", ufl_files)
@pytest.mark.parametrize("scalar_type", ["float64", "float32", "complex128", "complex64"])
def test_demo(file, scalar_type):
    """Test a demo."""
    if file in [
        "MixedGradient",
        "TraceElement",  # HDiv Trace
        "MixedElasticity",  # VectorElement of BDM
        "RestrictedElement",
        "_TensorProductElement",
    ]:
        # Skip demos that use elements not yet implemented in Basix
        pytest.skip()

    if "complex" in scalar_type and file in [
        "BiharmonicHHJ",
        "BiharmonicRegge",
        "StabilisedStokes",
    ]:
        # Skip demos that are not implemented for complex scalars
        pytest.skip()
    elif "Complex" in file and scalar_type in ["float64", "float32"]:
        # Skip demos that are only implemented for complex scalars
        pytest.skip()

    opts = f"--scalar_type {scalar_type}"
    extra_flags = "-Wunused-variable -Werror -fPIC "
    assert os.system(f"cd {demo_dir} && ffcx {opts} {file}.py") == 0
    assert (
        os.system(
            f"cd {demo_dir} && "
            "CPATH=../ffcx/codegeneration/ "
            f"gcc -I/usr/include/python{sys.version_info.major}.{sys.version_info.minor} "
            f"{extra_flags} "
            f"-shared {file}.c -o {file}.so"
        )
        == 0
    )
