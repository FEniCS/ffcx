import os
import sys

import pytest

demo_dir = os.path.dirname(os.path.realpath(__file__))

ufl_files = []
for file in os.listdir(demo_dir):
    if file.endswith(".py") and not file == "test_demos.py":
        ufl_files.append(file[:-3])


@pytest.mark.parametrize("file", ufl_files)
@pytest.mark.parametrize("scalar", ["float", "double"])
def test_demo(file, scalar):
    if file in [
        "MixedGradient",
        "TraceElement",  # HDiv Trace
        "MixedElasticity",  # VectorElement of BDM
        "RestrictedElement",
        "_TensorProductElement",
    ]:
        # Skip demos that use elements not yet implemented in Basix
        pytest.skip()

    extra_flags = "-Wunused-variable -Werror -fPIC "
    assert os.system(f"cd {demo_dir} && ffcx --scalar_type={scalar} {file}.py") == 0
    assert (
        os.system(
            f"cd {demo_dir} && "
            "CPATH=../ffcx/codegeneration/ "
            f"gcc -I/usr/include/python{sys.version_info.major}.{sys.version_info.minor} {extra_flags}"
            f"-shared {file}.c -o {file}.so"
        )
        == 0
    )


@pytest.mark.parametrize("file", ufl_files)
@pytest.mark.parametrize("scalar", ['"float _Complex"', '"double _Complex"'])
def test_demo_complex(file, scalar):
    if file not in [
        "CellGeometry",
        "VectorPoisson",
        "Symmetry",
        "ExpressionInterpolation",
        "MetaData",
        "MassDG0",
        "MixedCoefficient",
        "MathFunctions",
    ]:
        # Skip demos that fail with complex mode
        pytest.skip()

    extra_flags = "-Wunused-variable -Werror -fPIC "
    assert os.system(f"cd {demo_dir} && ffcx --scalar_type={scalar} {file}.py") == 0
    assert (
        os.system(
            f"cd {demo_dir} && "
            "CPATH=../ffcx/codegeneration/ "
            f"gcc -I/usr/include/python{sys.version_info.major}.{sys.version_info.minor} {extra_flags}"
            f"-shared {file}.c -o {file}.so"
        )
        == 0
    )
