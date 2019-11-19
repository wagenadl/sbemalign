
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

z0 = 804
nz = 5

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

sv = ri.subvolume(z0, nz)
mpp = subvol_elastic(sv)
subtract_mont(z0, mpp)
subtract_rigid(z0, mpp)
mp5.assignK(mpp)
ap = mp5.allpoints(mpp)
idx = mp5.elasticindex(mpp)
Ax, bx = mp5.elasticmatrix(mpp, idx, ap, 0)
dx = scipy.sparse.linalg.spsolve(Ax.tocsr(), bx)
dx = mp5.elasticdeindex(idx, dx)
xx=ap[(3,1,0)][0]
yy=ap[(3,1,0)][1]                   
idx1 = np.logical_and((xx-3256)**2 < 50**2, (yy+721)**2 < 50**2)
print(xx[idx1])
print(yy[idx1])
print(dx[(3,1,0)][idx1])
dx1 = scipy.sparse.linalg.spsolve(Ax.tocsr(), bx)

for mp in mpp:
    if mp.r1==2 and mp.r2==3 and mp.m2==1:
        idx1 = np.logical_and((mp.xx2-3256)**2 < 50**2, (mp.yy2+721)**2 < 50**2)
        k = np.nonzero(idx1)[0][0]
        p=idx[(3,1,0,k)]
        print(mp.xx2[idx1], mp.yy2[idx1], mp.xx2[idx1]-mp.xx1[idx1], mp.m1)
        print(p, Ax[p,p], bx[p], dx1[p])
        print(p-1, Ax[p-1,p-1], bx[p-1], dx1[p-1])

mp0 = None        
for mp in mpp:
    if mp.r1==2 and mp.r2==2 and mp.m1==0 and mp.m2==0 and mp.s2==32:
        mp0 = mp
        
        idx1 = np.logical_and((mp.xx2-3256)**2 < 100**2, (mp.yy2+721)**2 < 100**2)
        k = np.nonzero(idx1)[0][0]
        p=idx[(2,0,32,k)]
        print(mp.xx2[idx1], mp.yy2[idx1], mp.xx2[idx1]-mp.xx1[idx1], mp.m1)
        print(p, Ax[p,p], bx[p], dx1[p])
        
        
