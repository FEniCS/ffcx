# -*- coding: utf-8 -*-
import os, sys, glob
from os.path     import join, isfile, isdir, getsize, sep, normpath, dirname, basename
from distutils   import sysconfig
from swig_config import configure_swig_env, get_status_output

# Make sure that we have a good scons-version
EnsureSConsVersion(0, 96)

# Create a SCons Environment based on the main os environment
env = Environment(ENV=os.environ)

# UFC version number
major = 1
minor = 1

# Build the commandline options for SCons:
if env["PLATFORM"].startswith("win"):
    default_prefix = r"c:\local"
elif env["PLATFORM"] == "darwin":
    default_prefix = join(sep,"sw")
else:
    default_prefix = join(sep,"usr","local")

if env["PLATFORM"].startswith("win"):
    default_python_dir = os.path.join("$prefix", "Lib", "site-packages")
else:
    default_python_dir = sysconfig.get_python_lib(prefix="$prefix",plat_specific=True)

options = [
    # configurable options for installation:
    PathOption("prefix", "Installation prefix", default_prefix, PathOption.PathAccept),
    PathOption("includeDir", "ufc.h installation directory",
               join("$prefix","include"), PathOption.PathAccept),
    PathOption("pythonModuleDir", "Python module installation directory", 
               default_python_dir, PathOption.PathAccept),
    BoolOption("enablePyUFC", "Compile and install the python extension module", "Yes"),
    PathOption("pythonExtDir", "Python extension module installation directory",
               default_python_dir, PathOption.PathAccept),
    PathOption("boostDir", "Specify path to Boost", None),
    PathOption("pkgConfDir", "Directory for installation of pkg-config files",
               join("$prefix","lib","pkgconfig"), PathOption.PathAccept),
    BoolOption("cleanOldUFC", "Clean any old installed UFC modules", "No"),
    BoolOption("cacheOptions", "Cache command-line options for later invocations", "Yes")]

# Set the options using any cached options
cache_file = "options.cache"
opts = Options(cache_file, args=ARGUMENTS.copy())
opts.AddOptions(*options)
opts.Update(env)
cache_options = env.get("cacheOptions", False)
if cache_options:
    del env["cacheOptions"] # Don't store this value
    opts.Save(cache_file, env)
    # Restore cacheOptions
    env["cacheOptions"] = cache_options

env.Help(opts.GenerateHelpText(env))

# Start building the message presented at the end of a simulation
end_message = "\n"

# Check for old ufc installation
old_ufc_modules = []
for p in sys.path:
    if isfile(join(p,"ufc","__init__.py")):
        old_ufc_modules.append(join(p,"ufc"))

# If not cleaning
if not env.GetOption("clean"):
    # Notify the user about that options from scons/options.cache are being used:
    if isfile(cache_file) and getsize(cache_file):
        print "Using options from options.cache"
    
    # Append end_message if old ufc module excists
    if old_ufc_modules:
        end_message +=  """
---------------------------------------------------------
*** Warning: Old ufc module

%s
  
still excists in installation path.
Try remove these files with:

    scons -c cleanOldUFC=Yes

Note that you may need to be root.
"""%("\n".join("    " + m for m in old_ufc_modules))

# Generate pkgconfig file
pkg_config_file = "ufc-%d.pc" % major
file = open(pkg_config_file, "w")
file.write("Name: UFC\n")
file.write("Version: %d.%d\n" % (major, minor))
file.write("Description: Unified Form-assembly Code\n")
file.write("Cflags: -I%s\n" % normpath(env.subst(env["includeDir"])))
file.close()

# Set up installation targets
ufc_basename = join("src", "ufc", "ufc")
env.Install(env["pkgConfDir"],File(pkg_config_file))
env.Install(env["includeDir"],File(ufc_basename+".h"))
env.Install(join(env["pythonModuleDir"],"ufc_utils"),
            [File(f) for f in glob.glob(join("src","utils","python","ufc","*.py") )])

targets = [env["pkgConfDir"],env["pythonModuleDir"],env["includeDir"]]

# If compiling the extension module
if env["enablePyUFC"]:
    swig_env, message = configure_swig_env(env.Clone())
    end_message += message
    if swig_env["enablePyUFC"]:
        ufc_wrap, ufc_py = swig_env.CXXFile(target=ufc_basename,
                                            source=[ufc_basename+".i"])
        ufc_so = swig_env.SharedLibrary(target=ufc_basename, source=ufc_wrap)[0]
        
        # A SCons bug workaround. Fixed in newer SCons versions
        ufc_py = File(join(dirname(ufc_basename),basename(str(ufc_py))))

        # Set up installation targets
        env.Install(join(env["includeDir"], "swig"),File(ufc_basename+".i"))
        env.Install(env["pythonExtDir"],[ufc_py,ufc_so])
        targets.append(env["pythonExtDir"])

# Set the alias for install
env.Alias("install", targets)

# Create installation target folders if they don't exists:
if 'install' in COMMAND_LINE_TARGETS:
    for target_dir in [env.subst(d) for d in targets]:
        if not os.path.isdir(target_dir):
            os.makedirs(target_dir)

# If the user are cleaning remove any old ufc python modules
#if env.GetOption("clean") and isdir(join(env.subst(env["pythonModuleDir"]),"ufc")):
if old_ufc_modules and env["cleanOldUFC"]:
    Clean(env["pythonModuleDir"],old_ufc_modules)

if not env.GetOption("clean"):
    import atexit
    if 'install' not in COMMAND_LINE_TARGETS:
        end_message += """
---------------------------------------------------------
If there were no errors, run

    scons install

to install UFC on your system. Note that you may need
to be root in order to install. To specify an alternative
installation directory, run

    scons install prefix=<path>

---------------------------------------------------------
"""
    def out():
        print end_message
    atexit.register(out)

