"Collecting all commands here for automatic inclusion in cmdline utils."

from cmd_load     import add_load_options, run_load

from cmd_printing import add_str_options, run_str
from cmd_printing import add_repr_options, run_repr
from cmd_printing import add_tree_options, run_tree

from cmd_graphviz import add_graphviz_options, run_graphviz
from cmd_latex    import add_latex_options, run_latex

from cmd_analyse  import add_analyse_options, run_analyse

from cmd_compile  import add_compile_options, run_compile

from cmd_compile_dolfin  import add_compile_dolfin_options, run_compile_dolfin

from cmd_compile_pdelab import add_compile_pdelab_options, run_compile_pdelab

from cmd_test_code_formatting import \
    add_test_code_formatting_options, run_test_code_formatting

def get_version():
    return '0.1'

def add_default_options(opts):
    opts.add_option('-v', '--verbose', action='store', default='')
    opts.add_option('-d', '--debug', action='store', default='')
    opts.add_option('-q', '--quiet', action='store', default='')
    add_skipping_options(opts)

def add_skipping_options(opts):
    opts.add_option('--skip-elements', action='store_true', default=False)
    opts.add_option('--skip-coefficients', action='store_true', default=False)
    opts.add_option('--skip-expressions', action='store_true', default=False)
    opts.add_option('--skip-forms', action='store_true', default=False)
