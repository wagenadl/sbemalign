#!/usr/bin/python3

# Finds subvolumes mucked up by trans pug in matchpointsq5
# We need to find all subvolumes in which _any_ trans involved m1!=m2

import aligndb

db = aligndb.DB()
ri = db.runinfo()

for z0 in range(4, 9504, 100):
    print(z0)
    sv = ri.subvolume(z0, 200)
    rr = db.vsel(f'select distinct r from solveq5elastic where z0={z0}')[0]
    if len(sv) != rr.size:
        print('Mismatch: ', sv, ' vs ', rr)
