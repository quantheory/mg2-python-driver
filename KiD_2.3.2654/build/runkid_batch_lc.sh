#!/bin/bash

#MUSB -A climalg            # Allocation to charge job against                  
#MUSB -q pdebug             # Queue: batch, debug, or killable                  
#MUSB -N KiD_test           # Job name                                          
#MUSB -l nodes=1            # Compute Nodes (16 cores/node on cab)
#MUSB -l walltime=00:30:00  # Max run time (hh:mm:ss)
#MUSB -o KiD.%j.out         # stdout file (%j is jobID)                         
#MUSB -e KiD.%j.err         # stderr file (%j is jobID)                         

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
