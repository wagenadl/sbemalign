#!/usr/bin/python3

nz = 200

IX = IY = 5
X = Y = 684


import aligndb
import factory
import numpy as np
import matchpointsq5 as mp5

db = aligndb.DB()
ri = db.runinfo()

def subvol_mont(sv):
    mpp = []
    N = len(sv)
    for n in range(N):
        rl = sv[n]
        if n>0:
            mpp += mp5.MatchPoints.alltrans(rl.r-1, rl.r)
        mpp += mp5.MatchPoints.allcross(rl.r, rl.s0, rl.s1)
    return mpp

z0 = 7704
sv = ri.subvolume(z0, nz)
mpp = subvol_mont(sv)
idx = mp5.index(mpp)
Ax, bx = mp5.matrix(mpp, idx, 0)
Ay, by = mp5.matrix(mpp, idx, 1)
xm = np.linalg.solve(Ax, bx)
ym = np.linalg.solve(Ay, by)
tns = mp5.tension(mpp, idx, xm, ym)


