
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
