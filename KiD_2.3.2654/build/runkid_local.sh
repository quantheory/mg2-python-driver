#!/bin/bash

source /usr/global/tools/dotkit/init.sh  
use netcdf-intel-4.1.3

echo 'running KiD tests'
for f in *.nml; do
    echo "running KiD test: $f"
    ./KiD_1D.exe $f ${f%.*} > ${f%.*}.log
done

echo 'tests complete'
