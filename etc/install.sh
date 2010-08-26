#!/bin/bash
#
# Script for SLtools installation.
#
# 1) Compile some libraries written in C/C++ (Soon, Fortran files should appear)
# 2) Copy the package to another place and set the PATH env variable
#
# Checking dependencies for c/python modules creation
#

version=%VERSION%

SWIG=$(which swig 2> /dev/null )
CC=$(which cc 2> /dev/null )

if [ -z "$SWIG" -o -z "$CC" ]
then
    echo
    echo " Dependencies check error!"
    echo " Check if you have GCC compilers and Swig installed on your system."
    echo " Remember you also have to have Python headers (devel package) installed."
    echo
    exit 1
fi

echo ""
echo "Compile some C/Python libraries..."

# Operating System
#
OS=$(uname -s)

# Compile each function on their respective directories..
#
cd lens/get_nfw_concentration
make clean &> /dev/null

if [ "$OS" == "Linux" ]
then
    echo " Compiling get_nfw_concentration for $OS"
    make module_lnx &> /dev/null || { echo "Failed."; exit 1; }
elif [ "$OS" == "Darwin" ]
then
    echo " Compiling get_nfw_concentration for $OS"
    make module_mac &> /dev/null || { echo "Failed."; exit 1; }
else
    echo " Failed."
    exit 1
fi
cd - &> /dev/null


cd lens/compute_nfw_lens_parameters
make clean &> /dev/null

if [ "$OS" == "Linux" ]
then
    echo " Compiling compute_nfw_lens_parameters for $OS"
    make module_lnx &> /dev/null || { echo "Failed."; exit 1; }
elif [ "$OS" == "Darwin" ]
then
    echo " Compiling compute_nfw_lens_parameters for $OS"
    make module_mac &> /dev/null || { echo "Failed."; exit 1; }
else
    echo "Failed."
    exit 1
fi
cd - &> /dev/null

# ----------------------------------------------------------------------
# We can do better and install SLtools into another place on the system
#
echo ""
echo "Where do you want to install sltools? [${HOME}/lib]   "
echo "(Type Ctrl+C to cancel the installation if you want to)"
read INSTALLDIR

# Set the default value (HOME) nothing was given
#
[ -z "$INSTALLDIR" ] && INSTALLDIR="${HOME}/lib"

# Check if INSTALLDIR exist, if not, create it
#
[ ! -d "$INSTALLDIR" ] && { mkdir -p $INSTALLDIR || exit 1; }

# Check if it can be written
#
if [ ! -w "$INSTALLDIR" ]
then
    echo ""
    echo "ERROR: No write permission on $INSTALLDIR. Finishing $0"
    exit 1
fi

# Create directories where files will be copied, and do it
#


INSTALLDIR=${INSTALLDIR}/sltools-v$version

mkdir -p ${INSTALLDIR}/sltools/bin || exit 1

echo 
echo "Copying files..."
cp -r * ${INSTALLDIR}/sltools
rm -rf ${INSTALLDIR}/sltools/install.sh
echo

echo "------------------------------------------------------------------------"
echo ""
echo "Update your login file (.login, .bashrc or equivalent) with the following lines:"
echo ""
echo "   export PYTHONPATH=${INSTALLDIR}/sltools-v$version:\$PYTHONPATH"
echo "   export PATH=${INSTALLDIR}/sltools-v$version/sltools/bin:\$PATH"
echo ""
echo "------------------------------------------------------------------------"
echo "Done."
echo ""
