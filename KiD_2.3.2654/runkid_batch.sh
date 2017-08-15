#!/bin/bash

echo 'running KiD tests'
for f in *.nml; do
    echo "running KiD test: $f"
    ./runkid.sh $f > ${f%.*}.log &
    sleep 5
done
wait
echo 'tests complete'
