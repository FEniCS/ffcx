
def add_test_code_formatting_options(opts):
    pass

def run_test_code_formatting(options, args):
    from c_format_test import test_code_formatting
    from cpp_compiler import test_cpp_compilation
    return test_code_formatting() or test_cpp_compilation()
