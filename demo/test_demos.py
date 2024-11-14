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
    if sys.platform.startswith("win32") and "complex" in scalar_type:
        # Skip complex demos on win32
        pytest.skip(reason="_Complex not supported on Windows")

    if file in [
        "MixedGradient",
        "TraceElement",  # HDiv Trace
        "MixedElasticity",  # VectorElement of BDM
        "RestrictedElement",
        "_TensorProductElement",
    ]:
        # Skip demos that use elements not yet implemented in Basix
        pytest.skip(reason="Element not yet implemented in Basix")

    if "complex" in scalar_type and file in [
        "BiharmonicHHJ",
        "BiharmonicRegge",
        "StabilisedStokes",
    ]:
        # Skip demos that are not implemented for complex scalars
        pytest.skip(reason="Not implemented for complex types")
    elif "Complex" in file and scalar_type in ["float64", "float32"]:
        # Skip demos that are only implemented for complex scalars
        pytest.skip(reason="Not implemented for real types")

    if sys.platform.startswith("win32"):
        opts = f"--scalar_type {scalar_type}"
        extra_flags = "/std:c17"
        assert os.system(f"cd {demo_dir} && ffcx {opts} {file}.py") == 0
        assert (
            os.system(
                f"cd {demo_dir} && " f'cl.exe /I "../ffcx/codegeneration" {extra_flags} /c {file}.c'
            )
        ) == 0
        assert (
            os.system(
                f"cd {demo_dir} && "
                f'clang-cl.exe /I "../ffcx/codegeneration" {extra_flags} /c {file}.c'
            )
        ) == 0
    else:
        cc = os.environ.get("CC", "cc")
        opts = f"--scalar_type {scalar_type}"
        extra_flags = (
            "-std=c17 -Wunused-variable -Werror -fPIC -Wno-error=implicit-function-declaration"
        )
        assert os.system(f"cd {demo_dir} && ffcx {opts} {file}.py") == 0
        assert (
            os.system(
                f"cd {demo_dir} && "
                f"{cc} -I../ffcx/codegeneration "
                f"{extra_flags} "
                f"-c {file}.c"
            )
            == 0
        )


@pytest.mark.parametrize("scalar_type", ["float64", "float32"])
def test_demo_nvrtc(scalar_type):
    """Test generated CUDA code with NVRTC."""
    try:
        from nvidia import cuda_nvrtc
    except ImportError:
        pytest.skip(reason="Must have NVRTC pip package installed to run test.")

    if sys.platform.startswith("win32"):
        pytest.skip(reason="NVRTC CUDA wrappers not currently supported for Windows.")

    files = [
        "Components",
        "FacetIntegrals",
        "HyperElasticity",
        "MathFunctions",
        "StabilisedStokes",
        "VectorPoisson",
    ]
    opts = f"--scalar_type {scalar_type} --cuda_nvrtc"
    nvrtc_dir = os.path.dirname(os.path.realpath(cuda_nvrtc.__file__))
    cc = os.environ.get("CC", "cc")
    extra_flags = (
        "-std=c17 -Wunused-variable -Werror -fPIC -Wno-error=implicit-function-declaration"
    )
    for file in files:
        assert os.system(f"cd {demo_dir} && ffcx {opts} {file}.py") == 0
        assert (
            os.system(
                f"cd {demo_dir} && "
                f"{cc} -I../ffcx/codegeneration "
                f"{extra_flags} "
                f"-c {file}.c"
            )
            == 0
        )

    cxx = os.environ.get("CXX", "c++")
    assert (
        os.system(
            f"cd {demo_dir} && "
            f"{cxx} -I../ffcx/codegeneration -I{nvrtc_dir}/include -L{nvrtc_dir}/lib "
            f" -Werror -o nvrtc_test nvrtc_test.cpp "
            f"{' '.join([file+'.o' for file in files])} -l:libnvrtc.so.12"
        )
        == 0
    )
    assert os.system(f"LD_LIBRARY_PATH=$LD_LIBRARY_PATH:{nvrtc_dir}/lib {demo_dir}/nvrtc_test") == 0
