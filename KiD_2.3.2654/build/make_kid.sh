#!/bin/bash
# ==============================================================================
# Script to build the Kinemaitc Driver (KiD) source code on LLNL LC machines.
#
# Inputs (Optional):
#     1) Number of threads to use in parallel build
#     2) Where to copy executable and launch scripts
#
# D.J. Gardner @ LLNL
# Dec 2016
# ==============================================================================

# Set up LC enviornment
# Note: Also need to update src/compiler_options.inc
# --------------------------------------------------------------------
source /usr/global/tools/dotkit/init.sh # Enable dotkit packages
use icc-16.0.210                        # intel compiler
use netcdf-intel-4.1.3                  # netcdf for output
use python                              # Load python for use in makefile 
PYTHONHOME=/usr/apps/python/bin

# print software loaded with dotkit
use

# compile source code
# --------------------------------------------------------------------
BUILD_ROOT=$PWD
KID_ROOT=${HOME}/Climate/Physics/KiD/KiD_2.3.2654

cd $KID_ROOT
make -j $1 CASE=1D all

RESULT=$?
if [ ! $RESULT -eq 0 ]; then
    echo
    echo ">>> BUILD FAILED <<<"
    exit $RESULT
else
    echo ">>> BUILD COMPLETE <<<"
fi

# copy to testing directory
# --------------------------------------------------------------------
if [ $# -gt 1 ]; then 
    TESTDIR=$2
    if [ ! -d $TESTDIR ]; then
        mkdir -p $TESTDIR
    fi

    echo "Copying executable to $TESTDIR"
    cp $KID_ROOT/bin/KiD_1D.exe $TESTDIR/.

    echo "Copying batch launch script to $TESTDIR"
    cp $BUILD_ROOT/runkid_batch.sh $TESTDIR/.

    echo "Copying local launch script to $TESTDIR"
    cp $BUILD_ROOT/runkid_local.sh $TESTDIR/.

    echo "Copying local pattern launch script to $TESTDIR"
    cp $BUILD_ROOT/runkid_pattern.sh $TESTDIR/.
fi