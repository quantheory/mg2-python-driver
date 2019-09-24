#!/usr/bin/env bash

blocksize=4096
num_files=3
start_column=32768
last_file="$(($num_files - 1))"

for i in $(eval echo {0..$last_file})
do
    start="$(( $start_column + $i * $blocksize ))"
    end="$(( $start_column + ($i + 1) * $blocksize - 1))"
    echo "time jacobian3.py $start $end | tee jacobian_log_$start-$end.txt"
    nohup time jacobian3.py $start $end > jacobian_log_$start-$end.txt &
done
