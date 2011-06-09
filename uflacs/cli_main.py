"""
This module contains a reusable main function
that can detect commands from an app specific
module. Call it from a script with an app module
containing the following:

# Required for entire app:
def get_version():
    return '1.0'

# Optional for entire app:
def add_default_options(opts):
    pass

# Required for each command:
def run_COMMAND(options, args):
    return

# Optional for each command:
def add_COMMAND_options(opts):
    pass
"""

__author__ = "Martin Sandve Alnes"
__copyright__ = __author__

import optparse

def get_app_module_commands(app_module):
    return [n[4:] for (n, f) in vars(app_module).iteritems()
            if callable(f) and n.startswith("run_")]

def main(args, app_module):
    # Default to asking for help
    if not args:
        args = ['--help']

    # Get list of valid commands from application 
    cmds = get_app_module_commands(app_module)

    # Split command and arguments
    if args[0] in cmds:
        cmd, args = (args[0], args[1:]) 
    else:
        cmd = None
        # Allow various versions of asking for help
        if any(h in args for h in ('help', '-h', '--help')):
            # Syntax 'help command', '--help command', '-h command'
            if len(args) == 2 and args[1] in cmds:
                cmd, args = args[1], ['--help']
            else:
                args = ['--help']
        else:
            args = ['--help']

    # If we have a command, get its option and run functions
    if cmd:
        add_cmd_options = getattr(app_module, 'add_%s_options' % cmd)
        run_cmd = getattr(app_module, 'run_%s' % cmd)
        assert callable(add_cmd_options)
        assert callable(run_cmd)
    else:
        add_cmd_options = None
        run_cmd = None

    # Setup option parser with custom usage and version
    usage = '%prog <command> [options]'
    usage += '\n\nValid commands: '
    format = "\n  %%-%ds - %%s" % (max(map(len, cmds)),)
    for c in sorted(cmds):
        a = getattr(app_module, 'run_%s' % c)
        d = a.__doc__ or "<No short description.>"
        usage += format % (c, d)
    version = app_module.get_version()
    parser = optparse.OptionParser(usage=usage, version=version)

    # Setup default options shared between commands
    app_module.add_default_options(parser)

    # Setup command specific options in a separate group
    if cmd:
        title = 'Options for command %s' % (cmd,)
        doc = add_cmd_options.__doc__ or "<No args description.>"
        opts = optparse.OptionGroup(parser, title, doc)
        add_cmd_options(opts)
        parser.add_option_group(opts)

    # Parse argument list (eventually prints help text)
    options, args = parser.parse_args(args)

    # Execute command function
    if cmd:
        return run_cmd(options, args)
    else:
        return -1
