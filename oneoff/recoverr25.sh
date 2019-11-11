#!/bin/bash

python3 ./slicealignq5.py
python3 ./relmontalignq5.py
python3 ./relmontattouchq5.py

python3 ./optimizeq5rel.py
python3 ./warpq5run.py

python3 ./interrunq5.py
python3 ./transrunmontq5.py
