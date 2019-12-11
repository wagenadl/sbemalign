#!/usr/bin/python3

# Finds subvolumes mucked up by trans pug in matchpointsq5
# We need to find all subvolumes in which _any_ trans involved m1!=m2

import aligndb

db = aligndb.DB()
ri = db.runinfo()

zz0 = []

for z0 in range(4, 9504, 100):
    sv = ri.subvolume(z0, 200)
    got = False
    for k in range(1, len(sv)):
        r1 = sv[k-1].r
        r2 = sv[k].r
        cnt = db.sel(f'''select count(*) from transrunmontq5
        where ((r={r1} and r2={r2}) or (r={r2} and r2={r1}))
        and m!=m2''')[0][0]
        if cnt>0:
            got = True
            break
    if got:
        print(z0)
        zz0.append(z0)


for z0 in zz0:
    print(z0)
    db.exe(f'delete from solveq5elastic where z0={z0}')
    db.exe(f'delete from renderq5elasticdone where z>={z0} and z<{z0+200}')
    db.exe(f'delete from renderq1done where z>={z0} and z<{z0+200}')    

        
