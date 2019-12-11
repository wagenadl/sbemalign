#!/usr/bin/python3

# Finds subvolumes mucked up by trans pug in matchpointsq5
# We need to find all subvolumes in which _any_ trans involved m1!=m2

import aligndb
import numpy as np
import sys
db = aligndb.DB()
ri = db.runinfo()

zz0 = []
for z0 in range(4, 9504, 100):
    print(z0)
    bad = False
    sv = ri.subvolume(z0, 200)
    r,m,s = db.vsel(f'select distinct r,m,s from solveq5elastic where z0={z0}')
    for vv in sv:
        M = ri.nmontages(vv.r)
        for vs in range(vv.s0, vv.s1):
            N = np.sum(np.logical_and(r==vv.r, s==vs))
            if N!=M:
                print(f'Mismatch at Z{z0} R{vv.r} S{vs}: M={M} vs N={N}')
                bad = True
    if bad:
        zz0.append(z0)
print(zz0)

sys.exit(1)

for z in zz0:
    db.exe(f'delete from solveq5elastic where z0={z}')
    db.exe(f'delete from renderq5elasticdone where z>={z} and z<{z+200}')
    db.exe(f'delete from renderq1done where z>={z} and z<{z+200}')
    
