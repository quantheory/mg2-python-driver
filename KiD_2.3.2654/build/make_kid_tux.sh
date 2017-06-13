#!/bin/bash
# ==============================================================================
# Script to build the Kinemaitc Driver (KiD) source code on LLNL LC machines.
#
# D.J. Gardner @ LLNL
# Dec 2016
#
# modified by C.J. Vogl @ LLNL
# May 2017
# ==============================================================================

# compile source code
# --------------------------------------------------------------------
BUILD_ROOT=$PWD
KID_ROOT=${HOME}/workspace/kid/KiD_2.3.2654

cd $KID_ROOT
make -j $1 CASE=1D COMPILER=gfortran NCPATH=${HOME}/local all

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

    echo "Copying launch script to $TESTDIR"
    cp $BUILD_ROOT/runkid_batch_tux.sh $TESTDIR/.
fi
