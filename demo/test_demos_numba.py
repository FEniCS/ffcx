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
def test_demo(file):
    if file in [
        "MixedGradient",
        "TraceElement",  # HDiv Trace
        "MixedElasticity",  # VectorElement of BDM
        "RestrictedElement",
        "HyperElasticity",  # tables too big
        "BiharmonicRegge",
        "ComplexPoisson",
        "Mini",
        "_TensorProductElement",
    ]:
        # Skip demos that use elements not yet implemented in Basix
        pytest.skip()

    assert os.system(f"cd {demo_dir} && ffcx -L numba {file}.py") == 0
    module = importlib.import_module(f"{file}_numba")

    names = [getattr(module, w[0]) for w in inspect.getmembers(module) if "tabulate_tensor" in w[0]]
    for scalar_tp in [np.float32, np.float64]:
        wrapper = module.wrapper(scalar_tp, scalar_tp)
        compiled_functions = [wrapper(f) for f in names]
        print(compiled_functions)

    # Some demos failing with complex numbers
    if file not in ["MathFunctions", "Conditional", "FunctionOperators"]:
        scalar_tp = np.complex128
        real_tp = np.float64
        wrapper = module.wrapper(scalar_tp, real_tp)
        compiled_functions = [wrapper(f) for f in names]
        print(compiled_functions)

    os.unlink(f"{demo_dir}/{file}_numba.py")
