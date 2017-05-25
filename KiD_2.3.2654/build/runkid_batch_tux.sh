#!/bin/bash

echo 'running KiD tests'
for f in *.nml; do
    echo "running KiD test: $f"
    ./KiD_1D.exe $f ${f%.*} > ${f%.*}.log &
    sleep 10
done
wait
echo 'tests complete'
