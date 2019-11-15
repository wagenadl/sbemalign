#!/usr/bin/python3

# This optimizes the relative positions of all the tiles in a subvolume.
# It does elastic deformation, translating the individual points.
# See E&R p. 1615

monttbl = 'solveq5mont'
rigidtbl = 'solveq5rigidtile'
outtbl = 'solveq5elastic'
nz = 200

IX = IY = 5
X = Y = 684

nthreads = 12

import aligndb
import factory
import numpy as np
import matchpointsq5 as mp5
import scipy.sparse
import scipy.sparse.linalg

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
    y float,
    dx float,
    dy float )''')
    # x, y are the x_k on p. 1615.
    # dx, dy are the ξ_k on p. 1615.
    
def subvol_elastic(sv):
    mpp = []
    first = True
    for rl in sv:
        if first:
            first = False
        else:
            for m in range(ri.nmontages(rl.r)):
                mpp += mp5.MatchPoints.trans(rl.r, m,
                                             thr=7,
                                             perslice=True)
        C = ri.ncolumns(rl.r)
        R = ri.nrows(rl.r)
        for c in range(C):
            for r in range(R-1):
                mpp += mp5.MatchPoints.cross(rl.r, r*C+c, (r+1)*C+c,
                                             rl.s0, rl.s1,
                                             thr=7,
                                             perslice=True)
        for c in range(C-1):
            for r in range(R):
                mpp += mp5.MatchPoints.cross(rl.r, r*C+c, r*C+c+1,
                                             rl.s0, rl.s1,
                                             thr=7,
                                             perslice=True)
        for m in range(ri.nmontages(rl.r)):
            mpi = mp5.MatchPoints.intra(rl.r, m, rl.s0, rl.s1, thr=7)
            mpe = mp5.MatchPoints.edge(rl.r, m, rl.s0, rl.s1, thr=7)
            mpp += mp5.combine(mpi, mpe)
    return mpp

def subtract_mont(z0, mpp):
    r,m,x,y = db.vsel(f'select r,m,x,y from {monttbl} where z0={z0}')
    for mp in mpp:
        mp.xx1 += x[np.logical_and(r==mp.r1, m==mp.m1)]
        mp.yy1 += y[np.logical_and(r==mp.r1, m==mp.m1)]
        mp.xx2 += x[np.logical_and(r==mp.r2, m==mp.m2)]
        mp.yy2 += y[np.logical_and(r==mp.r2, m==mp.m2)]
    return mpp

def subtract_rigid(z0, mpp):
    r,m,s,x,y = db.vsel(f'select r,m,s,x,y from {rigidtbl} where z0={z0}')
    for mp in mpp:
        rm1 = np.logical_and(r==mp.r1, m==mp.m1)
        rms1 = np.logical_and(rm1, s==mp.s1)
        rm2 = np.logical_and(r==mp.r2, m==mp.m2)
        rms2 = np.logical_and(rm2, s==mp.s2)
        mp.xx1 += x[rms1]
        mp.yy1 += y[rms1]
        mp.xx2 += x[rms2]
        mp.yy2 += y[rms2]
    return mpp

def optisub(z0):
    sv = ri.subvolume(z0, nz)
    mpp = subvol_elastic(sv)
    subtract_mont(z0, mpp)
    subtract_rigid(z0, mpp)
    mp5.assignK(mpp)
    ap = mp5.allpoints(mpp)
    idx = mp5.elasticindex(mpp)
    Ax, bx = mp5.elasticmatrix(mpp, idx, ap, 0)
    Ay, by = mp5.elasticmatrix(mpp, idx, ap, 1)
    dx = scipy.sparse.linalg.spsolve(Ax.tocsr(), bx)
    dy = scipy.sparse.linalg.spsolve(Ay.tocsr(), by)
    dx = mp5.elasticdeindex(idx, dx)
    dy = mp5.elasticdeindex(idx, dy)
    return ap, dx, dy

def perhapsoptisub(z0):
    cnt = db.sel(f'select count(*) from {rigidtbl} where z0={z0}')
    if cnt[0][0]==0:
        print(f'Skipping {z0}+{nz} - not yet done at rigid level')
        return
    cnt = db.sel(f'select count(*) from {outtbl} where z0={z0}')
    if cnt[0][0]>0:
        return
    print(f'Working on {z0}+{nz}')
    ap, dx, dy = optisub(z0)
    with db.db:
        with db.db.cursor() as c:
            for rms,xy in ap.items():
                r,m,s = rms
                dx1 = dx[rms]
                dy1 = dy[rms]
                for k in range(len(xy[0])):
                    c.execute(f'''insert into {outtbl}
                    ( z0, r, m, s, x, y, dx, dy )
                    values
                    ( {z0}, {r}, {m}, {s},
                    {xy[0][k]}, {xy[1][k]}, {dx1[k]}, {dy1[k]} )''')
        
createtable()

fac = factory.Factory(nthreads)
R = ri.nruns()
Z = ri.z0(R) + ri.nslices(R)

for z0 in range(4, Z, nz//2):
    fac.request(perhapsoptisub, z0)
    
fac.shutdown()

