#!/usr/bin/python3

# This optimizes the relative positions of all the montages in a subvolume
# It does not optimize the positions of individual points.
# See E&R p. 1611

crosstbl = 'slicealignq5'
#intratbl = 'relmontalignq5'
#edgetbl  = 'relmontattouchq5'
#intertbl = 'interrunq5'
transtbl = 'transrunmontq5'

outtbl = 'stichq5mont'
IX = IY = 5
X = Y = 684
nz = 200

nthreads = 1

import aligndb
import time
import sys
import traceback

import swiftir
import pyqplot as qp
import numpy as np

db = aligndb.DB()
ri = db.runinfo()

def droptable():
    db.nofail(f'drop table {outtbl}')

def createtable():
    db.exe(f'''create table if not exists {outtbl} (
    z0 integer,
    r integer,
    m integer,
    x float,
    y float )''')

class MatchPoints:
    def __init__(self):
        self.r1 = None
        self.m1 = None
        self.s1 = None
        self.r2 = None
        self.m2 = None
        self.s2 = None
        self.xx1 = None
        self.yy1 = None
        self.xx2 = None
        self.yy2 = None
        def x(self, ax):
            # ax=0 for x or 1 for y
            # This is "x" on E&R p. 1611
            if ax:
                return self.yy1
            else:
                return self.xx1
        def xp(self, ax):
            # ax=0 for x or 1 for y
            # This is "x'" on E&R p. 1611
            if ax:
                return self.yy2
            else:
                return self.xx2

    def cross(r, m1, m2, s0=0, s1=None, thr=20, perslice=False):
        # Returns a single MatchPoints with combined data for all slices,
        # unless PERSLICE is True, in which case a list of MatchPoints
        # for individual slices is returned.
        swhere = f' and s>={s0}'
        if s1 is not None:
            swhere += f' and s<{s1}'
        (s, x1,y1, x2,y2, snr) = db.vsel(f'''select
        s,
        (ix1+0.5)*{X}+dx0/2-dx/2-dxb/2-dxc/2,
        (iy1+0.5)*{Y}+dy0/2-dy/2-dyb/2-dyc/2,
        (ix2+0.5)*{X}-dx0/2+dx/2+dxb/2+dxc/2,
        (iy2+0.5)*{Y}-dy0/2+dy/2+dyb/2+dyc/2,
        snrc from {crosstbl}
        where r={r} and m1={m1} and m2={m2} {swhere} and snrc>={thr}''')
        if perslice:
            mpp = []
            for s1 in np.unique(s):
                mp = MatchPoints()
                mp.r1 = mp.r2 = r
                mp.m1 = m1
                mp.m2 = m2
                mp.s1 = mp.s2 = s1
                mp.xx1 = x1[s==s1]
                mp.yy1 = y1[s==s1]
                mp.xx2 = x2[s==s1]
                mp.yy2 = y2[s==s1]
                mpp.append(mp)
            return mpp
        else:
            mp = MatchPoints()
            mp.r1 = mp.r2 = r
            mp.m1 = m1
            mp.m2 = m2
            mp.s1 = mp.s2 = None
            mp.xx1 = x1
            mp.yy1 = y1
            mp.xx2 = x2
            mp.yy2 = y2
            return mp
        
    def trans(r2, m2, thr=20, perslice=False):
        # Implicitly, r1=r2-1, s2=0, s1=S(r1)-1.
        # We return a list of MatchPoints, one for each existing m1.
        # If perslice is True, s1 and s2 are stored in the MatchPoints,
        # otherwise, None.
        (m1, x2,y2, x1,y1, snr) = db.vsel(f'''select
        m2,
        (ix+0.5)*{X}-dx/2-dxb/2,
        (iy+0.5)*{Y}-dy/2-dyb/2,
        x+dx/2+dxb/2,
        y+dy/2+dyb/2,
        snrb
        from {transtbl}
        where r={r2} and m={m2} and snrb>={thr}''')
        mm1 = np.unique(m1)
        mpp = []
        for m in mm1:
            mp = MatchPoints()
            mp.r1 = r2-1
            mp.m1 = m
            mp.r2 = r2
            mp.m2 = m2
            if perslice:
                mp.s1 = ri.nslices(r2-1) - 1
                mp.s2 = 0
            else:
                mp.s1 = mp.s2 = None
            mp.xx1 = x1[m1==m]
            mp.yy1 = y1[m1==m]
            mp.xx2 = x2[m1==m]
            mp.yy2 = y2[m1==m]
            mpp.append(mp)
        return mpp
    
    def subvol_mont(sv, crossthr, transthr):
        mpp = []
        first = True
        for rl in mpp:
            if first:
                first = False
            else:
                for m in range(ri.nmontages(rl.r)):
                    mpp += MatchPoints.trans(rl.r, m, transthr)
                    C = ri.ncolumns(rl.r)
                    R = ri.nrows(rl.r)
            for c in range(C):
                for r in range(R-1):
                    mpp.append(MatchPoints.cross(rl.r, r*C+c,
                                                 rl.r, (r+1)*C+c,
                                                 rl.s0, rl.s1,
                                                 crossthr))
            for c in range(C-1):
                for r in range(R):
                    mpp.append(MatchPoints.cross(rl.r, r*C+c,
                                                 rl.r, r*C+c+1,
                                                 rl.s0, rl.s1,
                                                 crossthr))
        return mpp

def index(mpp):
    # Given a list of MatchPoints objects, construct an index for matrixing
    k = 0
    idx = {} # Map from (r,m,s) to k
    for mp in mpp:
        k1 = (mp.r1, mp.m1, mp.s1)
        k2 = (mp.r2, mp.m2, mp.s2)
        if k1 not in idx:
            idx[k1] = k
            k += 1
        if k2 not in idx:
            idx[k2] = k
            k += 1
    return idx
    
def matrix(mpp, idx, q):
    EPSILON = 1e-6
    K = length(idx)
    A = np.eye(K) * EPSILON
    b = np.zeros(K)
    for mp in mpp:
        k = idx[(mp.r1, mp.m1, mp.s1)]
        kp = idx[(mp.r2, mp.m2, mp.s2)]
        w = len(mp.xx1) # Could be changed, of course
        Dx = np.mean(mp.xp(q) - mp.x(q))
        A[k,k] += w
        A[kp,kp] += w
        A[k,kp] -= w
        A[kp,k] -= w
        b[k] += Dx
        b[kp] -= Dx
    return A, b

def optisub(sv):
    mpp = MatchPoints.subvol_mont(sv, 20, 20)
    idx = index(mpp)
    Ax, bx = matrix(mpp, idx, 0)
    Ay, by = matrix(mpp, idx, 1)
    xm = np.linalg.solve(Ax, bx)
    ym = np.linalg.solve(Ay, by)

def perhapsoptisub(z0):
    cnt = db.sel(f'select count(*) from {outtbl} where z0={z0}')
    if cnt[0][0]==0:
        optisub(z0)
        
'''
createtable()

fac = factory.Factory(nthreads)
R = ri.nruns()
Z = ri.z0(R) + ri.nslices(R)

for z0 in range(4, Z, nz//2):
    fac.request(perhapsoptisub, z0)
    
fac.shutdown()
'''