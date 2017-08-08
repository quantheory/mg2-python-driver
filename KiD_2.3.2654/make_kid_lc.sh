#!/bin/bash
# ==============================================================================
# Script to build the Kinematic Driver (KiD) source code on LLNL LC machines.
#
# D.J. Gardner @ LLNL
# C.J. Vogl @ LLNL
# Jul 2016
# ==============================================================================

# Set up LC environnment
# --------------------------------------------------------------------
if [[ $HOSTNAME == "cab"* ]]; then
  SYSTEM=cab
  source /usr/global/tools/dotkit/init.sh # Enable dotkit packages
  use icc-16.0.210                        # intel compiler
  use netcdf-intel-4.1.3                  # netcdf for output
  use python                              # Load python for use in makefile 
  COMPILER=ifort
  NCPATH=/usr/local/tools/netcdf-intel-4.1.3
else
  SYSTEM=moodyblues
  COMPILER=gfortran
  NCPATH=/usr
fi

# compile source code
# --------------------------------------------------------------------
BUILD=build_$SYSTEM

make CASE=1D \
  COMPILER=$COMPILER \
  NCPATH=$NCPATH \
  EXECDIR=$BUILD/bin \
  OBJDIR=$BUILD/obj \
  all

RESULT=$?
if [ ! $RESULT -eq 0 ]; then
    echo
    echo ">>> BUILD FAILED <<<"
    exit $RESULT
else
    echo ">>> BUILD COMPLETE <<<"
fi
