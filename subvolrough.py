#!/usr/bin/python3

import aligndb
import time
import sys
import traceback

import numpy as np
import scipy.sparse
import scipy.sparse.linalg

db = aligndb.DB()
ri = db.runinfo()

X = Y = 684 # This is correct for our Q5 work
MAXX = X*5  # This is correct for our Q5 work
IX = IY = 5

crosstbl = 'slicealignq5'
intratbl = 'relmontalignq5'
edgetbl  = 'relmontattouchq5'
transtbl = 'transrunmontq5'

def roughtransdelta(r1, m1, r2, m2, thresh):
    '''Returns dx, dy, w.
    A pixel at (x,y) in (r1,m1) matches a pixel at (x-dx,y-dy) in (r2,m2).
    W is the number of matchpoints found.
    It is required that r2 = r1 + 1.
    If no matchpoints exist, returns 0, 0, 0.
    '''
    if r1 != r2-1:
        raise ValueError(f'roughtransdelta must have R{r2} != R{r1} + 1.')
    sdx,sdy,n = db.sel(f'''select 
    sum(x+dx+dxb-{X}*(ix+.5)), 
    sum(y+dy+dyb-{Y}*(iy+.5)), 
    count(*)
    from {transtbl}
    where r={r2} and m={m2} and m2={m1} and snrb>={thresh}''')[0]
    if n==0:
        return 0, 0, 0
    else:
        return sdx/n, sdy/n, n

def roughcrossdelta(r, m1, m2, thresh, s0=0, s1=None):
    '''Returns dx, dy, w.
    A pixel at (x,y) in (r,m1) matches a pixel at (x-dx,y-dy) in (r,m2).
    W is the number of matchpoints found.
    If no matchpoints exist, returns 0, 0, 0.
    '''
    if m2<m1:
        dx,dy,n = roughcrossdelta(r, m2, m1, thresh, s0, s1)
        return -dx, -dy, n
    if s1 is None:
        swhere = f's>={s0}'
    else:
        swhere = f's>=s0 and s<{s1}'
    sdx,sdy,n = db.sel(f'''select 
    sum(-dx-dxb-dxc-{X}*(ix2-ix1)+dx0), 
    sum(-dy-dyb-dyc-{Y}*(iy2-iy1)+dy0), 
    count(*)
    from {crosstbl}
    where r={r} and m1={m1} and m2={m2} and snrc>={thresh} and {swhere}''')[0]
    if n==0:
        return 0, 0, 0
    else:
        return sdx/n, sdy/n, n

class AllRoughDeltas:
    def __init__(self, ri, subvol, transthresh=20, crossthresh=20):
        self.cross = {} # keys are (r,m1,m2), values are (dx,dy,w)
        self.trans = {} # keys are (r1,m1,r2,m2), values are (dx,dy,w)
        N = len(subvol)
        for n in range(N):
            r = subvol[n].r
            M = ri.nmontages(r)
            C = ri.ncolumns(r)
            R = ri.nrows(r)
            for c in range(C-1):
                for ro in range(R):
                    m = c + C*ro
                    self.cross[(r,m,m+1)] = roughcrossdelta(r,m,m+1, crossthresh)
            for c in range(C):
                for ro in range(R-1):
                    m = c + C*ro
                    self.cross[(r,m,m+C)] = roughcrossdelta(r,m,m+C, crossthresh)
        for n in range(1,N):
            r1 = subvol[n-1].r
            r2 = subvol[n].r
            M2 = ri.nmontages(r2)
            for m2 in range(M2):
                mm1, = db.vsel(f'''select distinct m2 from {transtbl} 
                where r={r2} and m={m2}''')
                for m1 in mm1:
                    self.trans[(r1,m1,r2,m2)] = roughtransdelta(r1,m1,r2,m2,
                                                                transthresh)

def roughindex(ri, subvol):
    # Returns a mapping from (r,m) to 0..K-1
    idx = {}
    k=0
    N = len(subvol)
    for n in range(N):
        r = subvol[n].r
        M = ri.nmontages(r)
        for m in range(M):
            idx[(r,m)] = k
            k += 1
    return idx

class RoughEqns:
    def __init__(self, deltas, idx, coord):
        if coord=='x':
            self.which = 0
        elif coord=='y':
            self.which = 1
        else:
            raise ValueError(f'coord must be "x" or "y", not {coord}')

        self.deltas = deltas
        self.idx = idx
        self.K = len(self.idx)
        # We are setting up a matrix equation M x = b.
        # where x_k are the positions of the various montages.
        # We are minimizing
        #   χ² = Σ_rm ε (x_rm)²
        #        + Σ_rmr'm' α_rmr'm' (x_rm - x_r'm' - Δx_rmr'm')²,
        # where ε is a small number for numerical stability,
        # α is the weight from roughtransdelta or roughcrossdelta,
        # and Δx is the displacement from those functions.
        # (Note that I am writing "x", but "y" is analogous.)
        # There are K equations from ∂χ²/∂x_rm = 0. The terms of those
        # equations are filled by the "fill" functions.
        # Summarizing k = (r,m), I can write:
        #  ∂χ²/∂x_k = 2 ε x_k
        #             + 4 Σ_n α_k,n (x_k - x_n - Δx_k,n).
        # (I used the fact that α is symmetric and Δx antisymmetric.)
        # Thus A_kl = 2ε δ_kl + 4 (Σ_n α_k,n) δ_kl - 4 α_k,l
        # and b_k = 4 Σ_n α_k,n Δ_k,n.
        self.m = np.zeros((self.K, self.K))
        self.b = np.zeros(self.K)
        self.fillstable()
        self.fillcross()
        self.filltrans()

    def fillstable(self):
        # Implement the terms from Σ_rm ε (x_rm)².
        # If I summarize r,m as k, that part gives terms
        # ∂χ²/∂x_k = δ_kl ε x_l.
        eps = 1e-6
        for k in range(self.K):
            self.m[k,k] += eps

    def fillcross(self):
        for k,v in self.deltas.cross.items():
            r,m1,m2 = k
            dx,dy,w = v
            k1 = self.idx[(r,m1)]
            k2 = self.idx[(r,m2)]
            self.m[k1,k1] += w
            self.m[k2,k2] += w
            self.m[k1,k2] -= w
            self.m[k2,k1] -= w
            self.b[k1] -= w*v[self.which]
            self.b[k2] += w*v[self.which]

    def filltrans(self):
        for k,v in self.deltas.trans.items():
            r1,m1,r2,m2 = k
            dx,dy,w = v
            k1 = self.idx[(r1,m1)]
            k2 = self.idx[(r2,m2)]
            self.m[k1,k1] += w
            self.m[k2,k2] += w
            self.m[k1,k2] -= w
            self.m[k2,k1] -= w
            self.b[k1] -= w*v[self.which]
            self.b[k2] += w*v[self.which]

    def solve(self):
        x = np.linalg.solve(self.m, self.b)
        self.x = {}
        for rm,k in self.idx.items():
            self.x[rm] = x[k]
        return self.x
        
if __name__=='__main__':
    import rawimage
    import pyqplot as qp

    # Test roughtransdelta
    dx,dy,w = roughtransdelta(1,0, 2,0, 10)
    img1 = rawimage.partialq5img(1,0,ri.nslices(1)-1, 2, 2)
    img2 = rawimage.partialq5img(2,0,0, 2, 2)
    qp.figure('/tmp/s1', 8, 4)
    qp.subplot(1,2,1)
    qp.imsc(img1, xx=np.arange(X), yy=-np.arange(Y))
    qp.marker('+')
    qp.pen('r', 2)
    qp.mark(X/2+dx/2, -(Y/2+dy/2))
    qp.shrink(1,1)
    qp.subplot(1,2,2)
    qp.imsc(img2, xx=np.arange(X), yy=-np.arange(Y))
    qp.marker('+')
    qp.pen('r', 2)
    qp.mark(X/2-dx/2, -(Y/2-dy/2))
    qp.shrink(1,1)

    # Test roughcrossdelta
    dx,dy,w = roughcrossdelta(1,0,1, 10)
    img1 = rawimage.partialq5img(1,0,400, 2, 4)
    img2 = rawimage.partialq5img(1,1,400, 2, 0)
    dy -= 4*Y
    qp.figure('/tmp/s2', 8, 4)
    qp.subplot(1,2,1)
    qp.imsc(img1, xx=np.arange(X), yy=-np.arange(Y))
    qp.marker('+')
    qp.pen('r', 2)
    qp.mark(X/2+dx/2, -(Y/2+dy/2))
    qp.shrink(1,1)
    qp.subplot(1,2,2)
    qp.imsc(img2, xx=np.arange(X), yy=-np.arange(Y))
    qp.marker('+')
    qp.pen('r', 2)
    qp.mark(X/2-dx/2, -(Y/2-dy/2))
    qp.shrink(1,1)

    # Test all rough
    sv = ri.subvolume(700,200)
    ard = AllRoughDeltas(ri, sv)
    ridx = roughindex(ri, sv)
    xeqns = RoughEqns(ard, ridx, 'x')
    yeqns = RoughEqns(ard, ridx, 'y')
    xsol = xeqns.solve()
    ysol = yeqns.solve()

    # Just one example of cross
    r = 3
    m1 = 2
    m2 = 3
    ix1=4
    iy1=2
    ix2=0
    iy2=2
    s=50
    img1 = rawimage.partialq5img(r,m1,s,ix1,iy1)
    img2 = rawimage.partialq5img(r,m2,s,ix2,iy2)
    dx = xsol[(r,m1)] - xsol[(r,m2)] + X*(ix1-ix2)
    dy = ysol[(r,m1)] - ysol[(r,m2)] + Y*(iy1-iy2)
    qp.figure('/tmp/s3', 8, 4)
    qp.subplot(1,2,1)
    qp.imsc(img1, xx=np.arange(X), yy=-np.arange(Y))
    qp.marker('+')
    qp.pen('r', 2)
    qp.mark(X/2-dx/2, -(Y/2-dy/2))
    qp.shrink(1,1)
    qp.subplot(1,2,2)
    qp.imsc(img2, xx=np.arange(X), yy=-np.arange(Y))
    qp.marker('+')
    qp.pen('r', 2)
    qp.mark(X/2+dx/2, -(Y/2+dy/2))
    qp.shrink(1,1)
    
    # Just one example of trans
    r1 = 2
    m1 = 1
    r2 = 3
    m2 = 3
    # Following {ix,iy} are just to make sure the marks are inside the
    # subimages. Other choices would have been fine too.
    ix1=3
    iy1=2
    ix2=1
    iy2=4
    img1 = rawimage.partialq5img(r1,m1,ri.nslices(r1)-1,ix1,iy1)
    img2 = rawimage.partialq5img(r2,m2,0,ix2,iy2)
    dx = xsol[(r1,m1)] - xsol[(r2,m2)] + X*(ix1-ix2)
    dy = ysol[(r1,m1)] - ysol[(r2,m2)] + Y*(iy1-iy2)
    qp.figure('/tmp/s3', 8, 4)
    qp.subplot(1,2,1)
    qp.imsc(img1, xx=np.arange(X), yy=-np.arange(Y))
    qp.marker('+')
    qp.pen('r', 2)
    qp.mark(X/2-dx/2, -(Y/2-dy/2))
    qp.shrink(1,1)
    qp.subplot(1,2,2)
    qp.imsc(img2, xx=np.arange(X), yy=-np.arange(Y))
    qp.marker('+')
    qp.pen('r', 2)
    qp.mark(X/2+dx/2, -(Y/2+dy/2))
    qp.shrink(1,1)
