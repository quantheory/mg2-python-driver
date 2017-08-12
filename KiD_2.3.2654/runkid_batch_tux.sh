#!/bin/bash

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/local/netcdf-intel-4.1.3/lib
source /usr/apps/intel/15.5.223/composer_xe_2015.5.223/bin/iccvars.sh intel64
source /usr/apps/intel/15.5.223/composer_xe_2015.5.223/bin/ifortvars.sh intel64

echo 'running KiD tests'
for f in *.nml; do
    echo "running KiD test: $f"
    ./KiD_1D.exe $f ${f%.*} > ${f%.*}.log &
    sleep 10
done
wait
echo 'tests complete'
