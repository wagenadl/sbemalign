#!/usr/bin/python3

# Core facilities for optimizign the relative positions of all the
# montages in a subvolume at q5.
# See E&R p. 1611ff.

crosstbl = 'slicealignq5'
intratbl = 'relmontalignq5'
edgetbl  = 'relmontattouchq5'
intertbl = 'interrunq5'
transtbl = 'transrunmontq5'

IX = IY = 5
X = Y = 684

import aligndb
import numpy as np

db = aligndb.DB()
ri = db.runinfo()

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

    def combine(self, other):
        if self.r1 != other.r1 or self.r2 != other.r2:
            raise Exception('Cannot combine: run mismatch')
        if self.m1 != other.m1 or self.m2 != other.m2:
            raise Exception('Cannot combine: montage mismatch')
        if self.s1 != other.s1 or self.s2 != other.s2:
            raise Exception('Cannot combine: slice mismatch')
        self.xx1 = np.hstack((self.xx1, other.xx1))
        self.yy1 = np.hstack((self.yy1, other.yy1))
        self.xx2 = np.hstack((self.xx2, other.xx2))
        self.yy2 = np.hstack((self.yy2, other.yy2))

    def __repr__(self):
        n = f'-- {len(self.xx1)} pts'
        if self.r1==self.r2:
            if self.m1==self.m2:
                return f'R{self.r1} M{self.m1} S{self.s1}:{self.s2} {n}'
            else:
                if self.s1 is None:
                    return f'R{self.r1} M{self.m1}:{self.m2} {n}'
        if self.s1 is None:
            return f'R{self.r1}.M{self.m1} : R{self.r2}.M{self.m2} {n}'
        s = f'R{self.r1}.M{self.m1}.S{self.s1}'
        return s + f' : R{self.r2}.M{self.m2}.S{self.s2} {n}'

    def cross(r, m1, m2, s0=0, s1=None, thr=20, perslice=False):
        # Returns a single MatchPoints with combined data for all slices,
        # unless PERSLICE is True, in which case a list of MatchPoints
        # for individual slices is returned.
        swhere = f' and s>={s0}'
        if s1 is not None:
            swhere += f' and s<{s1}'
        #print(r, m1, m2, s0, s1, thr, perslice, swhere)
        try:
            (s, x1,y1, x2,y2, snr) = db.vsel(f'''select
            s,
            (ix1+0.5)*{X}+dx0/2-dx/2-dxb/2-dxc/2,
            (iy1+0.5)*{Y}+dy0/2-dy/2-dyb/2-dyc/2,
            (ix2+0.5)*{X}-dx0/2+dx/2+dxb/2+dxc/2,
            (iy2+0.5)*{Y}-dy0/2+dy/2+dyb/2+dyc/2,
            snrc from {crosstbl}
            where r={r} and m1={m1} and m2={m2} {swhere} and snrc>={thr}''')
        except Exception:
            print(f'Cross failed >= {thr} at R{r} M{m1}:{m2} S{s0}-{s1}')
            raise
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
        try:
            (m1, x2,y2, x1,y1, snr) = db.vsel(f'''select
            m2,
            (ix+0.5)*{X}-dx/2-dxb/2,
            (iy+0.5)*{Y}-dy/2-dyb/2,
            x+dx/2+dxb/2,
            y+dy/2+dyb/2,
            snrb
            from {transtbl}
            where r={r2} and m={m2} and snrb>={thr}''')
        except Exception:
            print(f'Trans failed to find any points >= {thr} at R{r2} M{m2}')
            raise
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
    
    def intra(r, m, s0, s1, thr):
        # Returns a list of MatchPoints with data from the intra table
        # for each of the slice pairs in [s0, s1).
        (s, x2,y2, x1,y1, snr) = db.vsel(f'''select
        s,
        (ix+0.5)*{X}-dx/2-dxb/2,
        (iy+0.5)*{Y}-dy/2-dyb/2,
        (ix+0.5)*{X}+dx/2+dxb/2,
        (iy+0.5)*{Y}+dy/2+dyb/2,
        snrb
        from {intratbl}
        where r={r} and m={m} and s>{s0} and s<{s1} and snrb>={thr}''')
        mpp = []
        for s2 in range(s0+1, s1):
            mp = MatchPoints()
            mp.r1 = mp.r2 = r
            mp.m1 = mp.m2 = m
            mp.s1 = s2 - 1
            mp.s2 = s2
            mp.xx1 = x1[s==s2]
            mp.yy1 = y1[s==s2]
            mp.xx2 = x2[s==s2]
            mp.yy2 = y2[s==s2]
            mpp.append(mp)
        return mpp

    def edge(r, m, s0, s1, thr):
        # Returns a list of MatchPoints with data from the edge table
        # for each of the slice pairs in [s0, s1).
        (s, x2,y2, x1,y1, snr) = db.vsel(f'''select
        s,
        (ix+0.5)*{X}+x-dx/2-dxb/2,
        (iy+0.5)*{Y}+y-dy/2-dyb/2,
        (ix+0.5)*{X}+x+dx/2+dxb/2,
        (iy+0.5)*{Y}+y+dy/2+dyb/2,
        snrb
        from {edgetbl}
        where r={r} and m={m} and s>{s0} and s<{s1} and snrb>={thr}''')
        mpp = []
        for s2 in range(s0+1, s1):
            mp = MatchPoints()
            mp.r1 = mp.r2 = r
            mp.m1 = mp.m2 = m
            mp.s1 = s2 - 1
            mp.s2 = s2
            mp.xx1 = x1[s==s2]
            mp.yy1 = y1[s==s2]
            mp.xx2 = x2[s==s2]
            mp.yy2 = y2[s==s2]
            mpp.append(mp)
        return mpp

def combine(mpi, mpe):
    mpr = []
    if len(mpi) != len(mpe):
        raise Exception('Combine needs equal length lists')
    for k in range(len(mpi)):
        mp = mpi[k]
        mp.combine(mpe[k])
        mpr.append(mp)
    return mpr
    
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

def deindex(idx, xx):
    res = {}
    for k,v in idx.items():
        res[k] = xx[v]
    return res
    
def matrix(mpp, idx, q):
    EPSILON = 1e-6
    K = len(idx)
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
        b[k] += w*Dx
        b[kp] -= w*Dx
    return A, b
