#!/bin/bash
# ==============================================================================
# Script to build the Kinemaitc Driver (KiD) source code on LLNL LC machines.
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

use

# Working directory and case name
# --------------------------------------------------------------------
WRKDIR=/p/lscratchd/${USER}/KiD
CASENAME='warm00_test'

# compile source code
# --------------------------------------------------------------------
BUILD_ROOT=$PWD
KID_ROOT=${HOME}/Climate/Physics/KiD/KiD_2.3.2654

cd $KID_ROOT
make clean
make -j $1 CASE=1D all

RESULT=$?
if [ $RESULT -eq 0 ]; then
    echo ">>> Successfully Built KiD <<<"
else
    echo ">>> Failed to Build KiD <<<"
    exit $RESULT
fi

# setup tests
# --------------------------------------------------------------------
RUNDIR=${WRKDIR}/${CASENAME}
if [ ! -d $RUNDIR ]; then
    mkdir -p $RUNDIR
fi

echo "Copying executable to $RUNDIR"
cp $KID_ROOT/bin/KiD_1D.exe $RUNDIR/.

echo "Copying launch script to $RUNDIR"
cp $BUILD_ROOT/runkid_batch.sh $RUNDIR/.

# clean up build
# ------------------------------------------------------------------------------
make clean