#!/bin/bash

#MSUB -A climalg            # Allocation to charge job against                  
#MSUB -q pdebug             # Queue: batch, debug, or killable                  
#MSUB -N KiD_test           # Job name                                          
#MSUB -l nodes=1            # Compute Nodes (16 cores/node on cab)
#MSUB -l walltime=00:30:00  # Max run time (hh:mm:ss)
#MSUB -o KiD.%j.out         # stdout file (%j is jobID)                         
#MSUB -e KiD.%j.err         # stderr file (%j is jobID)                         

source /usr/global/tools/dotkit/init.sh  
use netcdf-intel-4.1.3

echo 'running KiD tests'
for f in *.nml; do
    echo "running KiD test: $f"
    ./KiD_1D.exe $f ${f%.*} > ${f%.*}.log &
    sleep 10
done
wait
echo 'tests complete'
