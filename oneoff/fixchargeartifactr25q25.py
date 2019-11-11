#!/usr/bin/python3

import produceq25
import aligndb
import factory

db = aligndb.DB()
ri = db.runinfo()

r = 25
S = ri.nslices(r)
M = ri.nmontages(r)

fac = factory.Factory(12)
for m in range(M):
    for s in range(141, S):
        fac.request(produceq25.produceq25, r, m, s)
fac.shutdown()

        
