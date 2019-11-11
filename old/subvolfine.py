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

# Terminology:
# INTRA refers to alignment within a montage (adjacent s, same m)
# CROSS refers to alignment between montages (same s, adjacent m)
# EDGE refers to alignment within a montage at locations that are also
# used for cross alignment (adjacent s, same m, referencing adjacent m)
# TRANS refers to alignment between montages of different runs

class Delta:
    def __init__(self, SHP):
        # Data is organized as (S, NY, NX) for S slices and NYxNX overlap points
        self.ix = np.zeros(SHP, dtype=int) # subtile index
        self.iy = np.zeros(SHP, dtype=int) 
        self.xx = np.zeros(SHP) # point location in full slice
        self.yy = np.zeros(SHP)
        self.dx = np.zeros(SHP) # displacement of that point wrt model*
        self.dy = np.zeros(SHP)
        self.snr = np.zeros(SHP) # SNRC of the point
        # *: i.e., half of the full inter-montage displacement in the case
        #    of cross deltas, or the full displacement wrt the slice above
        #    in the case of intra-montage deltas
    def shape(self):
        return self.ix.shape
    def empty(self):
        return self.ix.size==0
    def __repr__(self):
        S, NY, NX = self.ix.shape
        return f'Delta[ {S} x {NY} x {NX} ]'

def crossdeltas(r, m1, m2, tbl, s0=0, s1=None):
    '''Finds the points of overlap between montages m and m_ in run r
    according to the named database table.
    Returns a pair of Delta structures.
    The first corresponds to points in m1, the second in m2.
    The first structure signifies that a pixel in (r,m1,s) at
    (xx,yy) corresponds to a pixel in (r,m2,s) at (xx+dx,yy+dy).
    '''
    S = ri.nslices(r)
    N = db.sel(f'''select count(*) from {tbl}
               where r={r} and m1={m1} and m2={m2} and s={s0}''')[0][0]
    NX = db.sel(f'''select count(distinct ix1) from {tbl}
                where r={r} and m1={m1} and m2={m2} and s={s0}''')[0][0]    
    NY = db.sel(f'''select count(distinct iy1) from {tbl}
                where r={r} and m1={m1} and m2={m2} and s={s0}''')[0][0]
    if NX*NY != N:
        raise Exception(f'Irregular contact matrix for R{r} M{m1}:{m2} in {tbl}')
    if s1 is None:
        s1 = S
    S = s1 - s0
    SHP = (S, NY, NX)
    res1 = Delta(SHP)
    res2 = Delta(SHP)
    if N==0:
        return (res1,res2)
    (s, ix1,iy1, ix2,iy2, dx0, dy0, dx,dy, snr) = db.vsel(f'''select
        s, ix1,iy1, ix2,iy2,
        dx0, dy0, dx+dxb+dxc,dy+dyb+dyc, snrc
        from {tbl}
        where r={r} and m1={m1} and m2={m2} and s>={s0} and s<{s1}
        order by s,iy1,ix1''')
    if len(s) != S * N:
        raise Exception(f'Mismatched point count for R{r} M{m1}:{m2} in {tbl}: {len(s)} rather than {S}*{N}={S*N}')

    x1 = ix1*X + X/2 + dx0/2 - dx/2
    x2 = ix2*X + X/2 - dx0/2 + dx/2
    y1 = iy1*Y + Y/2 + dy0/2 - dy/2
    y2 = iy2*Y + Y/2 - dy0/2 + dy/2
    # These positions are in the full tile
    dx = x2 - x1
    dy = y2 - y1

    res1.ix = np.reshape(ix1, SHP)
    res1.iy = np.reshape(iy1, SHP)
    res1.xx = np.reshape(x1, SHP)
    res1.yy = np.reshape(y1, SHP)
    res1.dx = np.reshape(dx/2, SHP)
    res1.dy = np.reshape(dy/2, SHP)
    res1.snr = np.reshape(snr, SHP)

    res2.ix = np.reshape(ix2, SHP)
    res2.iy = np.reshape(iy2, SHP)
    res2.xx = np.reshape(x2, SHP)    
    res2.yy = np.reshape(y2, SHP)    
    res2.dx = np.reshape(-dx/2, SHP)
    res2.dy = np.reshape(-dy/2, SHP)
    res2.snr = np.reshape(snr, SHP)
    return (res1, res2)

def _montagedeltas(r, where, tbl, name, xcol='x', ycol='y', s0=0, s1=None):    
    S = ri.nslices(r)
    if s1 is None:
        s1 = S
    S = s1 - s0
    N = db.sel(f'''select count(*) from {tbl}
               where {where} and s={s0}''')[0][0]
    NX = db.sel(f'''select count(distinct ix) from {tbl}
                where {where} and s={s0}''')[0][0]    
    NY = db.sel(f'''select count(distinct iy) from {tbl}
                where {where} and s={s0}''')[0][0]
    if NX*NY != N:
        raise Exception(f'Irregular contact matrix for {name} in {tbl}')
    SHP = (S, NY, NX)
    res = Delta(SHP)
    if N==0:
        return res
    (s, ix,iy, x,y, dx,dy, snr) = db.vsel(f'''select
        s, ix,iy, {xcol}, {ycol}, dx+dxb,dy+dyb, snrb
        from {tbl}
        where ({where}) and s>={s0} and s<{s1} order by s,iy,ix''')
    if len(s) != S*NY*NX:
        raise Exception(f'''Mismatched point count for {name} in {tbl}: 
        {len(s)} rather than {S}*{N}={S*N}''')
    res.xx = np.reshape(X*ix + x, SHP)
    res.yy = np.reshape(Y*iy + y, SHP)
    res.dx = np.reshape(dx, SHP)
    res.dy = np.reshape(dy, SHP)
    res.snr = np.reshape(snr, SHP)
    return res

def intradeltas(r, m, tbl, s0=0, s1=None):
    '''Finds the points where tiles in montage m in run r where compared
    to their (s+1 and s-1) according to the named database table.
    Returns a Delta structure.
    In the structure, (dx,dy) corresponds to the optimal shift of
    the corresponding point (xx,yy). That is, a pixel in (r,m,s) at
    (xx,yy) corresponds to a pixel at (xx+dx,yy+dy) in (r,m,s-1).'''
    return _montagedeltas(r, f'r={r} and m={m}', tbl, name=f'R{r} M{m}',
                          xcol=f'{X/2}-dx/2', ycol=f'{Y/2}-dy/2',
                          s0=s0, s1=s1)

def edgedeltas(r, m, m2, tbl, s0=0, s1=None):
    '''Finds the points where tiles in montage m in run r were compared
    to their neighbors (preceding and following slice) according to the
    named database table. Unlike INTRADELTAS, this function looks
    specifically for edge points, that is, points connected to a neighboring
    montage. (Called “touch points” elsewhere.)
    Returns a Delta structure.
    In the structure, (dx,dy) corresponds to the optimal shift of
    the corresponding point (xx,yy). That is, a pixel in (r,m,s) at
    (xx,yy) corresponds to a pixel at (xx+dx,yy+dy) in (r,m,s-1).
    The shape of the structure should match exactly to the first result
    of CROSSDELTAS(r, m, m2, tbl1) with the appropriate table.
    Note that the results here pertain to points in montage m, not in m2.'''
    return _montagedeltas(r, f'r={r} and m={m} and m2={m2}', tbl,
                          name=f'R{r} M{m}:{m2}',
                          xcol='x-dx/2', ycol='y-dy/2',
                          s0=s0, s1=s1)

class TransDelta(Delta):
    def __init__(self, SHP):
        # SHP must be 1xNYxNX
        super(TransDelta, self).__init__(SHP)
        self.m1 = np.zeros(SHP, dtype=int) # reference to montage in r-1.

def transdeltas(r, m, tbl):
    '''Finds the points of overlap between montage m in slice 0 of run r and 
    all montages in slice S-1 of run r-1. Returns a TransDelta structure.
    In that structure, xx,yy,ix,iy pertain to (r,m,0); dx,dy are the shift
    relative to (r-1,m1,S-1). That is, a pixel at (xx-dx, yy-dy) in (r,m0) 
    matches a pixel at (xx,yy) in (r-1,m1,S-1). '''
    ix,iy,m2,x,y,snr = db.vsel(f'''select 
    ix,iy,m2,x+dx+dxb,y+dy+dyb,snr
    from {tbl}
    where r={r} and m={m}
    order by iy,ix''')
    if len(ix) != IX*IY:
        raise Exception(f'Mismatch count in transdeltas: {len(ix)} != {IX}*{IY}')
    SHP = (IY,IX)
    res = TransDelta(SHP)
    ix = np.reshape(ix, SHP)
    iy = np.reshape(iy, SHP)
    m1 = np.reshape(m2, SHP)
    x = np.reshape(x, SHP)
    y = np.reshape(y, SHP)
    snr = np.reshape(snr, SHP)
    res.ix = ix
    res.iy = iy
    res.xx = X*(ix+.5)
    res.yy = Y*(iy+.5)
    res.dx = x - res.xx
    res.dy = y - res.yy
    res.m1 = m1
    res.snr = snr
    return res

class AllDeltas:
    '''AllDeltas are within one runlet. The top runlet has no trans.'''
    def __init__(self, r,
                 crosstbl=None, intratbl=None, edgetbl=None, transtbl=None,
                 s0=0, s1=None):
        self.r = r
        self.M = ri.nmontages(r)
        self.s0 = s0
        if s1 is None:
            s1 = ri.nslices(r)
        self.s1 = s1
        self.S = s1 - s0
        if crosstbl is not None:
            self.pullcross(crosstbl)
        if intratbl is not None:
            self.pullintra(intratbl)
        if edgetbl is not None:
            self.pulledge(edgetbl)
        if transtbl is not None:
            self.pulltrans(transtbl)

    def pullcross(self, tbl):
        self.cross = {}
        for m in range(self.M):
            for m_ in range(m, self.M):
                d1, d2 = crossdeltas(self.r, m, m_, tbl, self.s0, self.s1)
                self.cross[(m,m_)] = d1
                self.cross[(m_,m)] = d2

    def pullintra(self, tbl):
        self.intra = [] # internal points
        for m in range(self.M):
            self.intra.append(intradeltas(self.r, m, tbl, self.s0, self.s1))

    def pulledge(self, tbl):
        self.edge = []
        for m in range(self.M):
            delm = []
            for m_ in range(self.M):
                delm.append(edgedeltas(self.r, m, m_, tbl, self.s0, self.s1))
            self.edge.append(delm)

    def pulltrans(self, tbl):
        self.trans = []
        for m in range(self.M):
            self.trans.append(transdeltas(self.r, m, tbl))

'''
We will minimize an error function:

E = E_stable + E_intra + E_edge + E_cross

where

E_stable = ε Σ_k [x_k - x_k0]^2

(where x_k0 is the "naive" position of the point, which can be taken from the
intra, cross, or trans tables, depending on what kind of point "k" is.)

E_intra = β Σ_k [(x_k - x_k') - dxintra_k]^2

(where k=(r,s,m,ix,iy) and k'=(r,s-1,m,ix,iy) for s>0)

E_edge = β' Σ_k [(x_k - x_k') - dxedge_k]^2

(where k=(r,s,m,ii_m2) and k'=(r,s-1,m,ii_m2))

E_cross = γ Σ_k [(x_k - x_k') - dxcross_k]^2

(where k=(r,s,m,ii_m2) and k'=(r,s,m2,ii_m2))

E_trans = α Σ_k [(x_k - x_k') - dxtrans_k]^2

(where k=(r,s=0,m,ix,iy) and k'=(r-1,s=S-1,m',?,?)

I could also add elasticity terms

E_elast = ζ Σ_k [(x_k - x_k') - 684]^2 + ζ Σ_k [(x_k - x_k'')]^2

(where k=(r,s,m,ix,iy) and k'=(r,s,m,ix+1,iy) and k''=(r,s,m,ix,iy+1))

Solve by ∂E/∂x_p = 0. That yields P linear equations in the x_p's which
can be written in matrix form as A x = b.
A_pq is the coefficient of x_q in the p-th equation.
b_p is the constant in the p-th equation. (Mind the sign.)

For instance, ∂E_intra/∂x_k = 2 γ [x_k - x_k' - dxintra_k]
                              - 2 γ [x_k'' - x_k - dxintra_k'']
so there is a diagonal entry in A_kk: 4 γ,
an entry in b: 2 γ (dxintra_k - dxintra_k'').
There is also an off-diagonal term A_kk' = -2 γ
and an off-diagonal term A_kk'' = -2 γ.
Note that k' refers to the point above k (i.e., s-1) and k'' to the point 
below k (i.e., s+1).

Since x and y are completely independent, we'll solve them in separate steps.

First, we need to create an index so that all the various {x_k} can be put
in a single vector.

All of the above is slightly but possibly importantly wrong. That is because
none of the measurements are actually between equivalent points. In the
case of E_intra, this is subtle: I really did try to match up (r,m,s,ix,iy) to
(r,m,s-1,ix,iy), splitting the difference in dx. In the case of E_trans,
this is not at all subtle. I took a perfect point in (r,m,0,ix,iy) and found
whatever (r-1,m',S-1,x',y') to match to. That match point is not one of 
{(r-1,m',S-1,ix',iy')} at all. Meaning that it is not clear what diagonal
term to put in. I could yank at the nearest point, but perhaps I should yank
more gently on several nearby points?
For E_cross, the situation is also complicated. In my previous solution 
(optimizingrel.py), I decoupled the m and m' points, pulling each to half
the dx. But that is not really right. The correct approach, again, is to 
pull (r,m,s,ii_m') to some linear combination of actual points in the (r,m',s)
tile.
For E_edge, I think I can get away with pulling at the corresponding point
in s-1, because I split the differences. But this is not great; I should 
probably have stuck with predictable points rather than points based on 
wherever the slicealignq5 ended up measuring.
Overall, though, I think the only substantial issue is in E_trans. If I add
a strong enough elasticity term, I think it will be OK to pull on the nearest
point in E_trans.
To make all of this work, I think I'll have to pretend that my cross points
measure at regular points rather than at the real spots. That is not good.
Can I re-measure the cross points at the regular spots? No, because I only
have relative measures: There is not, in my approach, a "model" space to 
grab onto.
If I add elasticity terms not only for the intra points but for the whole
7x7 grid, I automatically regularize everything and my index becomes quite
simple.
'''

def fineindex(subvol):
    # Returns a mapping from (r,m,s,nx,ny) to 0..K-1.
    # Note that nx=1..IX correspond to the internal points; nx=0 and nx=IX+1
    # are the edge points. Same for ny. 
    idx = {}
    k = 0
    N = len(subvol)
    for n in range(N):
        r = subvol[n].r
        M = ri.nmontages(r)
        for s in range(subvol[n].s0, subvol[n].s1):
            for m in range(M):
                for nx in range(IX+2):
                    for ny in range(IY+2):
                        idx[r,m,s,nx,ny] = k
                        k = k + 1
    return idx

def FineEqns:
    def __init__(self, sv, deltas, idx, coord):
        '''SV must define a subvolume;
        DELTAS must be a list of ALLDELTAS, one per runlet in the subvolume;
        IDX must be from FINEINDEX;
        COORD must be either 'x' or 'y'.'''
        self.sv = sv
        self.N = len(sv)
        self.deltas = deltas
        self.idx = idx
        self.coord = coord
        if coord=='x':
            self.which = 0
        elif coord=='y':
            self.which = 1
        else:
            raise ValueError('coord must be "x" or "y"')
        
        self.K = len(self.idx)

        self.A = scipy.sparse.dok_matrix((K,K))
        self.b = np.zeros(self.K)
        self.fillstable(1e-6)
        self.fillelast(1e-2)
        #self.fillintra()
        #self.fillcross()
        #self.filledge()
        #self.filltrans()

    def fillstable(self, eps):
        for k in range(self.K):
            self.A[k,k] = eps
    def fillelast(self, zeta):
        '''Terms are ζ sum[k' neighbor of k] (x_k - x_k' - dx_natural_kk')^2.
        dE/dx_k = 2ζ [x_k - x_k' - dx_natural_kk']
                  -2ζ [x_k' - x_k - dx_natural_k'k]
'''
        for n in range(self.N):
            r = self.sv[n].r
            s0 = self.sv[n].s0
            s1 = self.sv[n].s1
            for m in range(ri.nmontages(r)):
                for s in range(s0, s1):
                    self.fillelast_internalh(r, m, s, zeta)
                    self.fillelast_internalv(r, m, s, zeta)
                    self.fillelast_edgeh(r, m, s, zeta)
                    self.fillelast_edgev(r, m, s, zeta)
    def fillelast_internalh(self, r, m, s, zeta):
        # Create internal horizontal edges
        for iy in range(IY):
            for ix in range(IX-1):
                k = self.idx[(r,m,s,ix+1,iy+1)]
                k_= self.idx[(r,m,s,ix+2,iy+1)]
                self.A[k,k] += zeta
                self.A[k_,k_] += zeta
                self.A[k,k_] -= zeta
                self.A[k_,k] -= zeta
                if self.which==0:
                    self.b[k] += X
                    self.b[k_] -= X
    def fillelast_internalv(self, r, m, s, zeta):
        # Create internal vertical edges
        for iy in range(IY-1):
            for ix in range(IX):
                k = self.idx[(r,m,s,ix+1,iy+1)]
                k_= self.idx[(r,m,s,ix+1,iy+2)]
                self.A[k,k] += zeta
                self.A[k_,k_] += zeta
                self.A[k,k_] -= zeta
                self.A[k_,k] -= zeta
                if self.which==1:
                    self.b[k] += Y
                    self.b[k_] -= Y
    def fillelast_edgeh(self, r, m, s, zeta):
        # Create horizontal edges between nx=0 and nx=1
        # and between nx=IX and nx=IX+1
        for ny in range(IY+2):
            for nx in [0, IX]:
                k = self.idx[(r,m,s,nx,ny)]
                k_= self.idx[(r,m,s,nx+1,ny)]
                self.A[k,k] += zeta
                self.A[k_,k_] += zeta
                self.A[k,k_] -= zeta
                self.A[k_,k] -= zeta
                if self.which==0:
                    self.b[k] += X*3/8
                    self.b[k_] -= X*3/8
    def fillelast_edgev(self, r, m, s, zeta):
        # Create vertical edges between ny=0 and ny=1
        # and between ny=IY and ny=IY+1
        for nx in range(IX+2):
            for ny in [0, IY]:
                k = self.idx[(r,m,s,nx,ny)]
                k_= self.idx[(r,m,s,nx,ny+1)]
                self.A[k,k] += zeta
                self.A[k_,k_] += zeta
                self.A[k,k_] -= zeta
                self.A[k_,k] -= zeta
                if self.which==1:
                    self.b[k] += Y*3/8
                    self.b[k_] -= Y*3/8
                           
