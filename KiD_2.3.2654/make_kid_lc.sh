#!/bin/bash
# ==============================================================================
# Script to build the Kinematic Driver (KiD) source code on LLNL LC machines.
#
# D.J. Gardner @ LLNL
# C.J. Vogl @ LLNL
# Sep 2017
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
elif [[ $HOSTNAME == "tux"* ]]; then
  SYSTEM=tux
  source /usr/apps/intel/15.5.223/composer_xe_2015.5.223/bin/iccvars.sh intel64
  source /usr/apps/intel/15.5.223/composer_xe_2015.5.223/bin/ifortvars.sh intel64
  COMPILER=ifort
  NCPATH=$HOME/local/netcdf-4.1.3_tux_intel_opt
  GPTLPATH=$HOME/local/gptl-5.5_tux_intel_opt
  MPIPATH=$HOME/local/mpich-3.2_tux_intel_opt
elif [[ $HOSTNAME == "MOODYBLUES" ]]; then
  SYSTEM=moodyblues
  COMPILER=gfortran
  NCPATH=/usr
else
  echo "Script not implement for this system"
  exit -1
fi

# set defines based on version number (v0 is default)
if [[ $1 == "v1" ]]; then
  export SED_UPDATECFL=True
elif [[ $1 == "v2" ]]; then
  export SED_UPDATECFL=True
  export SED_COMBINELAMBDA=True
elif [[ $1 == "v3" ]]; then
  export SED_UPDATECFL=True
  export SED_COMBINELAMBDA=True
  export SED_NONLINEAR=True
elif [[ $1 == "v4" ]]; then
  export SED_UPDATECFL=True
  export SED_COMPFLAG=True
elif [[ $1 == "v5" ]]; then
  export SED_UPDATECFL=True
  export SED_COMBINELAMBDA=True
  export SED_COMPFLAG=True
elif [[ $1 == "v6" ]]; then
  export SED_UPDATECFL=True
  export SED_COMBINELAMBDA=True
  export SED_COMPFLAG=True
  export SED_NONLINEAR=True
elif [[ $1 == "v7" ]]; then
  export SED_UPDATECFL=True
  export SED_COMBINELAMBDA=True
  export SED_USEWPA=True
elif [[ $1 == "v8" ]]; then
  export SED_UPATECFL=True
  export SED_COMBINELAMBDA=True
  export SED_NONLINEAR=True
  export SED_USEWPA=True
fi

# set build directory
BUILD=build_$SYSTEM
if [[ $2 == "new" ]]; then
  rm -rf $BUILD
fi
mkdir -p $BUILD

# compile source code
# --------------------------------------------------------------------
make CASE=1D \
  COMPILER=$COMPILER \
  NCPATH=$NCPATH \
  GPTLPATH=$GPTLPATH \
  MPIPATH=$MPIPATH \
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
