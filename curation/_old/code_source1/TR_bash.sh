#!/bin/bash
# script to launch code/change_TR.py in parallel
pwd

echo 'now starting parallel jobs'

N=1
(
while IFS= read -r path; do
    ((i=i%N)); ((i++==0)) && wait
    python code/change_TR.py "$path" &
done < "code/tr_edit.txt"
wait
)
