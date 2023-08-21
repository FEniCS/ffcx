import os
import numpy as np
import inspect
import importlib
import pytest

demo_dir = os.path.dirname(os.path.realpath(__file__))

ufl_files = []
for file in os.listdir(demo_dir):
    if file.endswith(".py") and "test" not in file:
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

    assert os.system(f"cd {demo_dir} && ffcx -L numba --scalar_type={scalar} {file}.py") == 0
    # assert (
    #     os.system(
    #         f"cd {demo_dir} && "
    #         "CPATH=../ffcx/codegeneration/ "
    #         f"gcc -I/usr/include/python{sys.version_info.major}.{sys.version_info.minor} {extra_flags}"
    #         f"-shared {file}.c -o {file}.so"
    #     )
    #     == 0
    # )
    module = importlib.import_module(f"{file}_numba")
    wrapper = module.wrapper(np.float64, np.float64)
    names = [getattr(module, w[0]) for w in inspect.getmembers(module) if "tabulate_tensor" in w[0]]
    compiled_functions = [wrapper(f) for f in names]
    print(compiled_functions)

# @pytest.mark.parametrize("file", ufl_files)
# @pytest.mark.parametrize("scalar", ['"float _Complex"', '"double _Complex"'])
# def test_demo_complex(file, scalar):
#     if file not in [
#         "CellGeometry",
#         "VectorPoisson",
#         "Symmetry",
#         "ExpressionInterpolation",
#         "MetaData",
#         "MassDG0",
#         "MixedCoefficient",
#         "MathFunctions",
#     ]:
#         # Skip demos that fail with complex mode
#         pytest.skip()

#     extra_flags = "-Wunused-variable -Werror -fPIC "
#     assert os.system(f"cd {demo_dir} && ffcx --scalar_type={scalar} {file}.py") == 0
#     assert (
#         os.system(
#             f"cd {demo_dir} && "
#             "CPATH=../ffcx/codegeneration/ "
#             f"gcc -I/usr/include/python{sys.version_info.major}.{sys.version_info.minor} {extra_flags}"
#             f"-shared {file}.c -o {file}.so"
#         )
#         == 0
#     )
