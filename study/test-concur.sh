#!/bin/zsh

for n in `seq 0 9`; do
    python3 ./test-concur.py $n &
done
