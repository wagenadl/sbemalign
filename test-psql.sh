#!/bin/zsh

for n in `seq 0 9`; do
    python3 ./test-psql.py $n &
done
