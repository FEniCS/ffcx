# Code generation format strings for UFC (Unified Form-assembly Code)
# This code is released into the public domain.
#
# The FEniCS Project (http://www.fenicsproject.org/) 2018
"""Code generation strings for an integral."""

declaration = """
extern ufcx_integral {factory_name};
"""

factory = """
// Code for integral {factory_name}

void tabulate_tensor_{factory_name}({scalar_type}* restrict A,
                                    const {scalar_type}* restrict w,
                                    const {scalar_type}* restrict c,
                                    const {geom_type}* restrict coordinate_dofs,
                                    const int* restrict entity_local_index,
                                    const uint8_t* restrict quadrature_permutation,
                                    void* custom_data)
{{
{tabulate_tensor}
}}

{enabled_coefficients_init}

ufcx_integral {factory_name} =
{{
  .enabled_coefficients = {enabled_coefficients},
  {tabulate_tensor_float32}
  {tabulate_tensor_float64}
  {tabulate_tensor_complex64}
  {tabulate_tensor_complex128}
  {tabulate_tensor_cuda_nvrtc}
  .needs_facet_permutations = {needs_facet_permutations},
  .coordinate_element_hash = {coordinate_element_hash},
  .domain = {domain},
}};

// End of code for integral {factory_name}
"""

cuda_wrapper = """

// Begin NVRTC CUDA wrapper for integral {factory_name}
// The wrapper is compiled with a standard C++ compiler, and is called at runtime to generate
// source code which is then compiled into a CUDA kernel at runtime via NVRTC.
void tabulate_tensor_cuda_nvrtc_{factory_name}(int* num_program_headers,
                                         const char*** program_headers,
                                         const char*** program_include_names,
                                         const char** out_program_src,
                                         const char** tabulate_tensor_function_name)
{{
  // The below typedefs are needed due to issues with including stdint.h in NVRTC source code
  const char* program_src = ""
    "#define alignas(x)\\n"
    "#define restrict __restrict__\\n"
    "\\n"
    "typedef unsigned char uint8_t;\\n"
    "typedef unsigned int uint32_t;\\n"
    "typedef double ufc_scalar_t;\\n"
    "\\n"
    "extern \\"C\\" __global__\\n"
    "void tabulate_tensor_{factory_name}({scalar_type}* restrict A,\\n"
    "                                    const {scalar_type}* restrict w,\\n"
    "                                    const {scalar_type}* restrict c,\\n"
    "                                    const {geom_type}* restrict coordinate_dofs,\\n"
    "                                    const int* restrict entity_local_index,\\n"
    "                                    const uint8_t* restrict quadrature_permutation\\n"
    "                                    )\\n"
    "{{\\n"
    "{tabulate_tensor_quoted}\\n"
    "}}";
  *num_program_headers = 0;
  *program_headers = NULL;
  *program_include_names = NULL;
  *out_program_src = program_src;
  *tabulate_tensor_function_name = "tabulate_tensor_{factory_name}";
}}

// End NVRTC CUDA wrapper for integral {factory_name}

"""


def get_factory(options):
    """Return the template string for constructing form integrals."""
    if options.get("cuda_nvrtc"):
        return cuda_wrapper + factory
    else:
        return factory
