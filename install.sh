#!/bin/bash

# ==========================================================================================
#
# Script for SLtools installation.
#
# Once SLtools has been downloaded and this package is unpacked, one just have to set up
# some files and SLtools should be ready to work.
#
# 1) Compile some libraries written in C/C++ (Soon, Fortran files should appear)
# 2) Copy the package to another place and set the PATH env variable
#
# ==========================================================================================


# Code version - VERSION is read from version.txt
. ./version.txt


# ------------------------------------------------------------------------------------------

# Checking dependencies for c/python modules creation
#
SWIG=$(which swig 2> /dev/null )
CC=$(which cc 2> /dev/null )

if [ -z "$SWIG" -o -z "$CC" ]
then
    echo
    echo " Dependencies check error!"
    echo " Check if you have GCC compiler and Swig installed on your system."
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
cd src/get_nfw_concentration
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


cd src/compute_nfw_lens_parameters
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

# Copy the C/Python compiled functions to their place on "tools" tree
#
cp src/get_nfw_concentration/module/*.py tools/simulations &> /dev/null
cp src/get_nfw_concentration/module/*.so tools/simulations &> /dev/null
cp src/get_nfw_concentration/module/*.conf tools/simulations &> /dev/null
cp src/compute_nfw_lens_parameters/module/*.py tools/lens &> /dev/null
cp src/compute_nfw_lens_parameters/module/*.so tools/lens &> /dev/null


# ----------------------------------------------------------------------
# We can do better and install SLtools into another place on the system
#
echo ""
echo "Where do you want to install SLtools? [${HOME}/lib]   "
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
DIR=${INSTALLDIR}/SLtools-${VERSION}
mkdir -p $DIR/sltools || exit 1
echo ""
echo "Copying files to $DIR/sltools ..."
cp -R * ${DIR}/sltools/.
mv ${DIR}/sltools/bin ${DIR}/sltools/.
echo ""


#[ -z "$PYTHONPATH" ] && PYTHONPATH=${DIR} || PYTHONPATH=${DIR}:${PYTHONPATH}
#export PYTHONPATH

echo "If you want SLtools to be readable in your PYTHON environment, use the command-line below."
echo "Update your login file (.login, .profile, ...) with the following line for permanent use."
echo ""
echo "   export PYTHONPATH=${DIR}:\$PYTHONPATH"
echo "   export PATH=${DIR}/bin:\$PATH"
echo ""
