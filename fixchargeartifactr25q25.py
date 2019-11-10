#!/usr/bin/python3

import produceq25.py
import aligndb

db = aligndb.DB()
ri = db.runinfo()

r = 25
S = ri.nslices(r)
M = ri.nmontages(r)
for m in range(M):
    for s in range(S):
        produceq25.produceq25(r, m, s)
        
