import setuptools
import sys

if sys.version_info < (3, 10):
    # Can be removed when pip editable user installs are fixed
    # https://github.com/pypa/pip/issues/7953
    import site
    site.ENABLE_USER_SITE = "--user" in sys.argv[1:]

setuptools.setup()
