
def add_test_code_formatting_options(opts):
    pass

def run_test_code_formatting(options, args):
    from uflacs.codeutils.cpp_format_test import test_cpp_formatting
    from uflacs.codeutils.cpp_compiler import test_cpp_compilation
    return test_cpp_formatting() or test_cpp_compilation()

