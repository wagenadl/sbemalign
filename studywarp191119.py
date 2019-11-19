#!/usr/bin/python3

import

import aligndb
import numpy as np

monttbl = 'solveq5mont'
rigidtbl = 'solveq5rigidtile'
elastbl = 'solveq5elastic'
globtbl = 'q5global'
bboxtbl = 'q5bbox'

X = Y = 684*5 # Full tile size in q5 space!

db = aligndb.DB()
ri = db.runinfo()

r = 21
s = 50

tiles = []
for m in ri.nmontages(r):
    tiles[m] = rawimage.
