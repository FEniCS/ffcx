
def add_example_options(opts):
    "See more examples at http://www.doughellmann.com/PyMOTW/optparse/"
    # Example string option:
    opts.add_option('-b', '--bing', action='store')

    # Example int option:
    opts.add_option('-c', '--count', action='store', type='int')

    # Example float option:
    opts.add_option('-d', '--delta', action='store', type='float')

    # Example bool option:
    opts.add_option('-a', '--all', action='store_true', default=False)

    # Example append option:
    opts.add_option('-o', '--object', action='append', dest='objects', default=[])

    # Example choice option:
    opts.add_option('-e', '--enum', type='choice', choices=('a', 'b', 'c'))

    # Example choice by option name:
    opts.add_option('--earth', action="store_const", const='earth', dest='element', default='earth')
    opts.add_option('--air',   action='store_const', const='air',   dest='element')
    opts.add_option('--water', action='store_const', const='water', dest='element')
    opts.add_option('--fire',  action='store_const', const='fire',  dest='element')

    # Example counting argument occurences:
    opts.add_option('-v', action="count", dest='verbosity', default=1)
    opts.add_option('-q', action='store_const', const=0, dest='verbosity')

