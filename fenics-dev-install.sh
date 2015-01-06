#!/usr/bin/env bash
#
# Developer script for configure + build + rebuild.
#
# Notes:
#
# - This script is what most developers use to build/rebuild this package.
# - This script is common to all distutils based FEniCS packages.
# - If this script is updated in one package, please propagate to the others!
#
# Environment variables:
#
# - $PROCS                    : controls number of processes to use for build
#                             : defaults to 6
# - $FENICS_PYTHON_EXECUTABLE : name of python executable
#                             : defaults to "python"
# - $FENICS_INSTALL_PREFIX    : path to FEniCS installation prefix
#                             : defaults to "${HOME}/opt/<branchname>"
#
# Note: Some of the code below may be redundant for either distutils or
# CMake based installations but it helps keeping the scripts up-to-date
# if the different scripts share as much code as possible.

# Exit on first error
set -e

# Get branch name
BRANCH=`(git symbolic-ref --short HEAD 2> /dev/null || git describe HEAD) | sed s:/:.:g`
echo "On branch '${BRANCH}'."

# Get installation prefix
: ${FENICS_INSTALL_PREFIX:="${HOME}/opt/fenics/${BRANCH}"}
echo "Installation prefix set to '${FENICS_INSTALL_PREFIX}'."

# Get Python executable and version
: ${FENICS_PYTHON_EXECUTABLE:=python}
FENICS_PYTHON_VERSION=$(${FENICS_PYTHON_EXECUTABLE} -c 'import sys; print(".".join(map(str, sys.version_info[:2])))')
echo "Python executable and version set to '${FENICS_PYTHON_EXECUTABLE} ${FENICS_PYTHON_VERSION}'."

# Build
${FENICS_PYTHON_EXECUTABLE} setup.py build

# Install
${FENICS_PYTHON_EXECUTABLE} setup.py install --prefix=${FENICS_INSTALL_PREFIX}

# Write config file
CONFIG_FILE="${FENICS_INSTALL_PREFIX}/fenics.conf"
rm -f ${CONFIG_FILE}
cat << EOF > ${CONFIG_FILE}
# FEniCS configuration file created by fenics-dev-install.sh on $(date)
export FENICS_INSTALL_PREFIX=${FENICS_INSTALL_PREFIX}
export FENICS_PYTHON_EXECUTABLE=${FENICS_PYTHON_EXECUTABLE}
export FENICS_PYTHON_VERSION=${FENICS_PYTHON_VERSION}

# Common Unix variables
export LD_LIBRARY_PATH=\${FENICS_INSTALL_PREFIX}/lib:\${LD_LIBRARY_PATH}
export PATH=\${FENICS_INSTALL_PREFIX}/bin:\${PATH}
export PKG_CONFIG_PATH=\${FENICS_INSTALL_PREFIX}/pkgconfig:\${PKG_CONFIG_PATH}
export PYTHONPATH=\${FENICS_INSTALL_PREFIX}/lib/python${FENICS_PYTHON_VERSION}/site-packages:\${PYTHONPATH}
export MANPATH=\${FENICS_INSTALL_PREFIX}/share/man:\${MANPATH}

# Cmake search path
export CMAKE_PREFIX_PATH=\${FENICS_INSTALL_PREFIX}:\${CMAKE_PREFIX_PATH}
EOF
if [ $(uname) = "Darwin" ]; then
    cat << EOF >> $CONFIG_FILE
# Mac specific path
export DYLD_FALLBACK_LIBRARY_PATH=\${FENICS_INSTALL_PREFIX}:\${DYLD_FALLBACK_LIBRARY_PATH}
EOF
fi

# Print information
echo
echo "- Installed branch '${BRANCH}' to ${FENICS_INSTALL_PREFIX}."
echo
echo "- Config file written to ${CONFIG_FILE}"
echo "  (source this file)."
echo
