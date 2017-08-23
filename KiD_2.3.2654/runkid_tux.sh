#!/bin/bash

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/local/netcdf-4.1.3_tux_intel_opt/lib
source /usr/apps/intel/15.5.223/composer_xe_2015.5.223/bin/iccvars.sh intel64
source /usr/apps/intel/15.5.223/composer_xe_2015.5.223/bin/ifortvars.sh intel64

rm -f ${1%.*}.*.txt
./KiD_1D.exe $1 ${1%.*}
