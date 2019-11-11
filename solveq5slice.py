#!/usr/bin/python3

# This optimizes the relative positions of all the tiles in a slice.
# It does not optimize the positions of individual points.
# See E&R p. 1611ff.
# Requires only slicealignq5 as input.

outtbl = 'solveq5slice'
crossthr = 20

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
    r integer,
    m integer,
    s integer,
    x float,
    y float )''')

def subvol_slice(r, s):
    mpp = []
    C = ri.ncolumns(r)
    R = ri.nrows(r)
    for c in range(C):
        for row in range(R-1):
            mpp += mp5.MatchPoints.cross(r, row*C+c, (row+1)*C+c,
                                         s, s+1,
                                         crossthr,
                                         True)
        for c in range(C-1):
            for row in range(R):
                mpp += mp5.MatchPoints.cross(r, row*C+c, row*C+c+1,
                                             s, s+1,
                                             crossthr,
                                             True)
    return mpp

def optislice(r, s):
    mpp = subvol_slice(r, s)
    idx = mp5.index(mpp)
    Ax, bx = mp5.matrix(mpp, idx, 0)
    Ay, by = mp5.matrix(mpp, idx, 1)
    xm = np.linalg.solve(Ax, bx)
    ym = np.linalg.solve(Ay, by)
    xm = mp5.deindex(idx, xm)
    ym = mp5.deindex(idx, ym)
    return xm, ym

def perhapsoptislice(r, s, force=False):
    if not force:
        cnt = db.sel(f'select count(*) from {outtbl} where r={r} and s={s}')
        if cnt[0][0]>0:
            return
    print(f'Working on R{r} S{s}')
    xm, ym = optislice(r, s)
    with db.db:
        with db.db.cursor() as c:
            c.execute(f'delete from {outtbl} where r={r} and s={s}')
            for k in xm:
                print(r, k, s, xm[k], ym[k])
                c.execute(f'''insert into {outtbl}
                ( r, m, s, x, y )
                values
                ( {r}, {k[1]}, {s}, {xm[k]}, {ym[k]} )''')
        
createtable()

fac = factory.Factory(nthreads)

for r0 in range(ri.nruns()):
    r = r0+1
    fac.request(perhapsoptislice, r, 0)
    S = ri.nslices(r)
    if S>1:
        fac.request(perhapsoptislice, r, S-1)
    
fac.shutdown()

