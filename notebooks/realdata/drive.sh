#!/bin/bash
cd /home/jrich/Desktop/fastQpick/notebooks/realdata
B=200
seq 1 $B | xargs -P 8 -I{} bash -c '
  s=$(printf "%04d" {})
  python run_one.py {} tmp/w{} counts/rep_${s}.txt 2>/dev/null && echo "done {}" >> progress.log
'
echo "ALL_DONE" >> progress.log
