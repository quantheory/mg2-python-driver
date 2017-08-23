#!/bin/bash

echo 'running KiD tests'
rm -f *-full.nml
for f in *.nml; do
  echo "running KiD test: $f"
  { time ./runkid.sh $f ; } &> ${f%.*}.log &
  sleep 5
done
wait
echo 'tests complete'
