
monttbl = 'solveq5mont'
rigidtbl = 'solveq5rigidtile'

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
    db.exe('create index if not exists sq5e_rms on solveq5elastic (r, m, s)')

def subvol_elastic(sv):
    mpp = []
    first = True
    for rl in sv:
        if first:
            first = False
        else:
            mpp += mp5.MatchPoints.alltrans(rl.r-1, rl.r,
                                             thr=7,
                                             perslice=True)
        mpp += mp5.MatchPoints.allcross(rl.r, rl.s0, rl.s1, thr=7, perslice=True)
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

nz = 3
z0 = ri.z0(20) + 68
z0a = 3404

sv = ri.subvolume(z0, nz)
mpp = subvol_elastic(sv)
subtract_mont(z0a, mpp)
subtract_rigid(z0a, mpp)
mp5.assignK(mpp)
ap = mp5.allpoints(mpp)
idx = mp5.elasticindex(mpp)
Ax, bx = mp5.elasticmatrix(mpp, idx, ap, 0)
Ay, by = mp5.elasticmatrix(mpp, idx, ap, 1)
dx = scipy.sparse.linalg.spsolve(Ax.tocsr(), bx)
dy = scipy.sparse.linalg.spsolve(Ay.tocsr(), by)
dx = mp5.elasticdeindex(idx, dx)
dy = mp5.elasticdeindex(idx, dy)

x,y=ap[(20,2,69)]
dx1=dx[(20,2,69)]
import pyqplot as qp
qp.figure('/tmp/s1')
qp.mark(y[x>1650], dx1[x>1650])
qp.yaxis()

qp.figure('/tmp/s2')
qp.marker('+')
qp.mark(x, y)

xx1, yy1, dx1, dy1 = db.vsel(f'''select x, y, dx, dy from solveq5elastic
   where r=20 and m=2 and s=69 and z0=3404''')
qp.figure('/tmp/s3')
qp.marker('+')
qp.mark(xx1,yy1)

