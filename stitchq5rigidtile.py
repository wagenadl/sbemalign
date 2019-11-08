#!/usr/bin/python3

# This optimizes the relative positions of all the tiles in a subvolume.
# It just does rigid translation.
# It does not optimize the positions of individual points.
# See E&R p. 1613

intbl = 'stitchq5mont'
outtbl = 'stitchq5rigidtile'
crossthr=20
transthr=15
transthr_despair=10
intrathr=20
edgethr=10
nz = 200

IX = IY = 5
X = Y = 684

nthreads = 12

import aligndb
import factory
import numpy as np
import matchpointsq5 as mp5

db = aligndb.DB()
ri = db.runinfo()

def droptable():
    db.nofail(f'drop table {outtbl}')

def createtable():
    db.exe(f'''create table if not exists {outtbl} (
    z0 integer,
    r integer,
    m integer,
    s integer,
    x float,
    y float )''')
    
def subvol_rigidtile(sv):
    mpp = []
    first = True
    for rl in sv:
        if first:
            first = False
        else:
            for m in range(ri.nmontages(rl.r)):
                try:
                    mpp += mp5.MatchPoints.trans(rl.r, m, transthr,
                                                 perslice=True)
                except:
                    mpp += mp5.MatchPoints.trans(rl.r, m, transthr_despair,
                                                 perslice=True)
        C = ri.ncolumns(rl.r)
        R = ri.nrows(rl.r)
        for c in range(C):
            for r in range(R-1):
                mpp += mp5.MatchPoints.cross(rl.r, r*C+c, (r+1)*C+c,
                                             rl.s0, rl.s1,
                                             crossthr,
                                             perslice=True)
        for c in range(C-1):
            for r in range(R):
                mpp += mp5.MatchPoints.cross(rl.r, r*C+c, r*C+c+1,
                                             rl.s0, rl.s1,
                                             crossthr,
                                             perslice=True)
        for m in range(ri.nmontages(rl.r)):
            mpi = mp5.MatchPoints.intra(rl.r, m, rl.s0, rl.s1, intrathr)
            mpe = mp5.MatchPoints.edge(rl.r, m, rl.s0, rl.s1, edgethr)
            mpp += mp5.combine(mpi, mpe)
    return mpp

def subtract_mont(z0, mpp):
    r,m,x,y = db.vsel(f'select r,m,x,y from {intbl} where z0={z0}')
    for mp in mpp:
        mp.xx1 += x[np.logical_and(r==mp.r1, m==mp.m1)]
        mp.yy1 += y[np.logical_and(r==mp.r1, m==mp.m1)]
        mp.xx2 += x[np.logical_and(r==mp.r2, m==mp.m2)]
        mp.yy2 += y[np.logical_and(r==mp.r2, m==mp.m2)]
    return mpp

def optisub(z0):
    sv = ri.subvolume(z0, nz)
    mpp = subvol_rigidtile(sv)
    subtract_mont(z0, mpp)
    idx = mp5.index(mpp)
    Ax, bx = mp5.matrix(mpp, idx, 0)
    Ay, by = mp5.matrix(mpp, idx, 1)
    xm = np.linalg.solve(Ax, bx)
    ym = np.linalg.solve(Ay, by)
    xm = mp5.deindex(idx, xm)
    ym = mp5.deindex(idx, ym)
    return xm, ym

def perhapsoptisub(z0):
    cnt = db.sel(f'select count(*) from {intbl} where z0={z0}')
    if cnt[0][0]==0:
        print(f'Skipping {z0}+{nz} - not yet done at mont level')
        return
    cnt = db.sel(f'select count(*) from {outtbl} where z0={z0}')
    if cnt[0][0]>0:
        return
    print(f'Working on {z0}+{nz}')
    xm, ym = optisub(z0)
    with db.db:
        with db.db.cursor() as c:
            for k in xm:
                c.execute(f'''insert into {outtbl}
                ( z0, r, m, s, x, y )
                values
                ( {z0}, {k[0]}, {k[1]}, {k[2]}, {xm[k]}, {ym[k]} )''')
        
createtable()

fac = factory.Factory(nthreads)
R = ri.nruns()
Z = ri.z0(R) + ri.nslices(R)

for z0 in range(4, Z, nz//2):
    fac.request(perhapsoptisub, z0)
    
fac.shutdown()

