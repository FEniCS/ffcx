import setuptools

# Can be removed when pip editable user installs are fixed
# https://github.com/pypa/pip/issues/7953
import site
import sys
site.ENABLE_USER_SITE = "--user" in sys.argv[1:]

setuptools.setup()
