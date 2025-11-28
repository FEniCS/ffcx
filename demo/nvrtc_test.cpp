#include "Components.h"
#include "FacetIntegrals.h"
#include "HyperElasticity.h"
#include "MathFunctions.h"
#include "StabilisedStokes.h"
#include "VectorPoisson.h"
#include "ufcx.h"
#include "nvrtc.h"
#include <iostream>
#include <stdexcept>
#include <sstream>
#include <string>
#include <vector>

void check_nvrtc_compilation(ufcx_form* form)
{
  // extract kernel
  ufcx_integral* integral = form->form_integrals[0];
  ufcx_tabulate_tensor_cuda_nvrtc* kernel = integral->tabulate_tensor_cuda_nvrtc;
  // call kernel to get CUDA-wrapped source code
  int num_program_headers;
  const char** program_headers;
  const char** program_include_names;
  const char* program_src;
  const char* tabulate_tensor_function_name;
  if (!kernel) {
    throw std::runtime_error("NVRTC wrapper function is NULL!");
  }
  (*kernel)(
    &num_program_headers, &program_headers,
    &program_include_names, &program_src,
    &tabulate_tensor_function_name);
  // compile CUDA-wrapped source code with NVRTC
  // with proper error checking

  nvrtcResult nvrtc_err;
  nvrtcProgram program;
  nvrtc_err = nvrtcCreateProgram(
    &program, program_src, tabulate_tensor_function_name,
    num_program_headers, program_headers,
    program_include_names);

 if (nvrtc_err != NVRTC_SUCCESS) {
    throw std::runtime_error(
      "nvrtcCreateProgram() failed with " +
      std::string(nvrtcGetErrorString(nvrtc_err)) + " "
      "at " + std::string(__FILE__) + ":" + std::to_string(__LINE__));
  }

  int num_compile_options = 0;
  const char** compile_options;
  // Compile the CUDA C++ program
  nvrtcResult nvrtc_compile_err = nvrtcCompileProgram(
    program, num_compile_options, compile_options);
  if (nvrtc_compile_err != NVRTC_SUCCESS) {
    // If the compiler failed, obtain the compiler log
    std::string program_log;
    size_t log_size;
    nvrtc_err = nvrtcGetProgramLogSize(program, &log_size);
    if (nvrtc_err != NVRTC_SUCCESS) {
      program_log = std::string(
        "nvrtcGetProgramLogSize() failed with " +
        std::string(nvrtcGetErrorString(nvrtc_err)) + " "
        "at " + std::string(__FILE__) + ":" + std::to_string(__LINE__));
    } else {
      program_log.resize(log_size);
      nvrtc_err = nvrtcGetProgramLog(
        program, const_cast<char*>(program_log.c_str()));
      if (nvrtc_err != NVRTC_SUCCESS) {
        program_log = std::string(
          "nvrtcGetProgramLog() failed with " +
          std::string(nvrtcGetErrorString(nvrtc_err))) + " "
          "at " + std::string(__FILE__) + ":" + std::to_string(__LINE__);
      }
      if (log_size > 0)
        program_log.resize(log_size-1);
    }
    nvrtcDestroyProgram(&program);

    std::stringstream ss;
    ss << "nvrtcCompileProgram() failed with "
       << nvrtcGetErrorString(nvrtc_compile_err) << "\n"
       << "CUDA C++ source code:\n"
       << std::string(60, '-') << "\n"
       << program_src
       << std::string(60, '-') << "\n"
       << "NVRTC compiler log:\n"
       << std::string(60, '-') << "\n"
       << program_log << "\n"
       << std::string(60, '-') << "\n";
    throw std::runtime_error(ss.str());
  }
}

int main()
{
  std::vector<ufcx_form*> forms = {
    form_Components_L,
    form_FacetIntegrals_a,
    form_HyperElasticity_a_F, form_HyperElasticity_a_J,
    form_MathFunctions_a,
    form_StabilisedStokes_a, form_StabilisedStokes_L,
    form_VectorPoisson_a, form_VectorPoisson_L  
  };
  
  for (ufcx_form* form : forms) check_nvrtc_compilation(form);
 
  return 0;
}

