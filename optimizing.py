#!/usr/bin/python3

import aligndb
import time
import sys
import traceback

import pyqplot as qp
import numpy as np
import scipy.sparse
import scipy.sparse.linalg

db = aligndb.DB()
ri = db.runinfo()

X = Y = 684 # This is correct for our Q5 work
MAXX = X*5  # This is correct for our Q5 work

# Terminology:
# INTRA refers to alignment within a montage (adjacent s, same m)
# CROSS refers to alignment between montages (same s, adjacent m)
# EDGE refers to alignment within a montage at locations that are also
# used for cross alignment (adjacent s, same m, referencing adjacent m)

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
        #    of cross deltas, or the full displacement wrt the stack average
        #    in the case of intraage deltas
    def shape(self):
        return self.ix.shape
    def empty(self):
        return self.ix.size==0
    def __repr__(self):
        S, NY, NX = self.ix.shape
        return f'Delta[ {S} x {NY} x {NX} ]'

def crossdeltas(r, m1, m2, tbl):
    '''Finds the points of overlap between montages m and m_ in run r
    according to the named database table.
    Returns a pair of Delta structures.
    The first corresponds to points in m1, the second in m2.
    In each structure, (dx,dy) corresponds to the optimal shift of
    the corresponding point (xx,yy). That is, a pixel in (r,m1,s) at
    (xx-dx,yy-dy)₁ corresponds to a pixel in (r,m2,s) at (xx-dx,yy-dy)₂.
    '''
    S = ri.nslices(r)
    N = db.sel(f'''select count(*) from {tbl}
               where r={r} and m1={m1} and m2={m2} and s=0''')[0][0]
    NX = db.sel(f'''select count(distinct ix1) from {tbl}
                where r={r} and m1={m1} and m2={m2} and s=0''')[0][0]    
    NY = db.sel(f'''select count(distinct iy1) from {tbl}
                where r={r} and m1={m1} and m2={m2} and s=0''')[0][0]
    if NX*NY != N:
        raise Exception(f'Irregular contact matrix for R{r} M{m1}:{m2} in {tbl}')
    SHP = (S, NY, NX)
    res1 = Delta(SHP)
    res2 = Delta(SHP)
    if N==0:
        return (res1,res2)
    (s, ix1,iy1, ix2,iy2, dx0, dy0, dx,dy, snr) = db.vsel(f'''select
        s, ix1,iy1, ix2,iy2,
        dx0, dy0, dx+dxb+dxc,dy+dyb+dyc, snrc
        from {tbl}
        where r={r} and m1={m1} and m2={m2} order by s,iy1,ix1''')
    if len(s) != S * N:
        raise Exception(f'Mismatched point count for R{r} M{m1}:{m2} in {tbl}')

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

def montageposition(deltas):
    '''Given a list of matrix of crossdeltas, determine optimal positions.
    This function handles one dimension at a time.
    DELTAS must be an MxM-matrix with nan entries where no data exists.
    Result is an M-vector.'''
    # We want to optimize E0 = sum_(defined pairs ij) (x_i - x_j - delta_ij)^2
    # For stability, we add to that E1 = epsilon sum_i x_i^2.
    # Note that
    # dE0/dx_i = 2 sum_(defined pairs with i) (x_i - x_j - delta_ij)
    # and that
    # dE1/dx_i = 2 epsilon x_i
    # Thus, we can solve this with a matrix eqn. of the form A x - b = 0,
    # where A_ii = 2 epsilon + 2 sum_j(pair ij exists)
    # A_ij = - 2 I(pair ij exists)
    # b_i = 2 sum_j delta_ij I(pair ij exists)
    # (I am confused about the factor two, but I know it is the same confusion
    # for A and b, so it doesn't matter except for the scale of E1.)
    M = deltas.shape[0]
    A = np.zeros((M,M)) + .001 # that's our epsilon
    b = np.zeros(M)
    for m in range(M):
        for m_ in range(M):
            if not np.isnan(deltas[m,m_]):
                A[m,m] += 1
                A[m,m_] -= 1
                b[m] += deltas[m,m_] - deltas[m_,m]
    x = np.linalg.solve(A, b)
    return x

def _montagedeltas(r, where, tbl, name, xcol='x', ycol='y'):    
    S = ri.nslices(r)
    N = db.sel(f'''select count(*) from {tbl}
               where {where} and s=0''')[0][0]
    NX = db.sel(f'''select count(distinct ix) from {tbl}
                where {where} and s=0''')[0][0]    
    NY = db.sel(f'''select count(distinct iy) from {tbl}
                where {where} and s=0''')[0][0]
    if NX*NY != N:
        raise Exception(f'Irregular contact matrix for {name} in {tbl}')
    SHP = (S, NY, NX)
    res = Delta(SHP)
    if N==0:
        return res
    (s, ix,iy, x,y, dx,dy, snr) = db.vsel(f'''select
        s, ix,iy, {xcol}, {ycol}, dx+dxb,dy+dyb, snrb
        from {tbl}
        where {where} order by s,iy,ix''')
    if len(s) != S*NY*NX:
        raise Exception(f'Mismatched point count for {name} in {tbl}')
    res.xx = np.reshape(X*ix + x, SHP)
    res.yy = np.reshape(Y*iy + y, SHP)
    res.dx = np.reshape(dx, SHP)
    res.dy = np.reshape(dy, SHP)
    res.snr = np.reshape(snr, SHP)
    return res

def intradeltas(r, m, tbl):
    '''Finds the points where tiles in montage m in run r where compared
    to their (s+1 and s-1) according to the named database table.
    Returns a Delta structure.
    In the structure, (dx,dy) corresponds to the optimal shift of
    the corresponding point (xx,yy). That is, a pixel in (r,m,s) at
    (xx-dx,yy-dy) corresponds to a pixel at (xx,yy) in montage-global
    coordinates.'''
    return _montagedeltas(r, f'r={r} and m={m}', tbl, name=f'R{r} M{m}',
                          xcol=f'{X/2}-dx', ycol=f'{Y/2}-dy')

def edgedeltas(r, m, m2, tbl):
    '''Finds the points where tiles in montage m in run r were compared
    to their neighbors (preceding and following slice) according to the
    named database table. Unlike INTRADELTAS, this function looks
    specifically for edge points, that is, points connected to a neighboring
    montage. (Called “touch points” elsewhere.)
    Returns a Delta structure.
    In the structure, (dx,dy) corresponds to the optimal shift of
    the corresponding point (xx,yy). That is, a pixel in (r,m,s) at
    (xx-dx,yy-dy) corresponds to a pixel at (xx,yy) in montage-global
    coordinates.
    The shape of the structure should match exactly to the first result
    of CROSSDELTAS(r, m, m2, tbl1) with the appropriate table.
    Note that the results here pertain to points in montage m, not in m2.'''
    return _montagedeltas(r, f'r={r} and m={m} and m2={m2}', tbl,
                          name=f'R{r} M{m}:{m2}',
                          xcol='x', ycol='y')

class AllDeltas:
    def __init__(self, r, crosstbl=None, intratbl=None, edgetbl=None):
        self.r = r
        self.M = ri.nmontages(r)
        self.S = ri.nslices(r)
        if crosstbl is not None:
            self.pullcross(crosstbl)
        if intratbl is not None:
            self.pullintra(intratbl)
        if edgetbl is not None:
            self.pulledge(edgetbl)

    def pullcross(self, tbl):
        self.cross = {}
        for m in range(self.M):
            for m_ in range(m, self.M):
                d1, d2 = crossdeltas(self.r, m, m_, tbl)
                self.cross[(m,m_)] = d1
                self.cross[(m_,m)] = d2

    def makemontpos(self):
        avgshiftx = np.zeros((self.M,self.M)) + np.nan
        avgshifty = np.zeros((self.M,self.M)) + np.nan
        for m in range(self.M):
            for m_ in range(self.M):
                crs = self.cross[(m,m_)]
                n = crs.dx.size
                if n>0:
                   avgshiftx[m,m_] = np.mean(crs.dx)
                   avgshifty[m,m_] = np.mean(crs.dy)
        posx = montageposition(avgshiftx)
        posy = montageposition(avgshifty)
        self.montpos = []
        for m in range(self.M):
            self.montpos.append((posx[m], posy[m]))

    def pullintra(self, tbl):
        self.intra = [] # internal points
        for m in range(self.M):
            self.intra.append(intradeltas(self.r, m, tbl))

    def pulledge(self, tbl):
        self.edge = []
        for m in range(self.M):
            delm = []
            for m_ in range(self.M):
                delm.append(edgedeltas(self.r, m, m_, tbl))
            self.edge.append(delm)

    '''
We will minimize an error function:

E = E_elast + E_intra + E_edge + E_cross

where

E_elast = α Σ_ksm [dx_ksm]^2

E_intra = β Σ_k,s>0,m [(dx_ksm - dx_k,s-1,m) - (dxintra_ksm - dxintra_k,s-1,m)]^2
E_edge = β' Σ_t,s>0,m [(dx_tsm - dx_t,s-1,m) - (dxedge_tsm - dxedge_t,s-1,m)]^2

E_cross = γ Σ_t,s,m,m' [dx_tsm - dxcross_ksmm']^2

Solve by dE/d dx_p = 0. That yields P linear equations in the dx_p's which
can be written in matrix form as A dx - b = 0.
A_pq is the coefficient of dx_q in the p-th equation.
b_p is the constant in the p-th equation. (Mind the sign.)

For instance, dE_cross/d dx_tsmm = 2gamma [dx_tsmm - deltaxstar_ksmm]
so there is a diagonal entry in A: 2gamma, and an entry in b: 2gamma deltaxstar.
Entries need to be added from the various terms in E.

Since x and y are completely independent, we'll solve them in separate steps.

First, we need to construct an index system to convert between k,m / t,m / t,m,m'
and p indices.
    '''

class Index:
    def __init__(self, ad,
                 crosssnrthr = 20,
                 intrasnrthr = 20,
                 edgesnrthr = 20):
        # Let's construct a linear indexing scheme for all our variables
        # We'll put the x coordinates first, followed by all the y coords.
        p = 0
        (self.intra, p) = self.makeintraindex(ad, p, intrasnrthr)
        (self.edge, p) = self.makeedgeindex(ad, p, crosssnrthr, edgesnrthr)
        self.P = p

    def makeedgeindex(self, ad, p, crossthr, edgethr):
        edge = []
        for m in range(ad.M):
            idxt = []
            for m_ in range(ad.M):
                mm_ = (m,m_)
                S,NY,NX = ad.edge[m][m_].shape()
                idx = np.zeros((S,NY,NX), dtype=int) - 1
                for ny in range(NY):
                    for nx in range(NX):
                        for s in range(S):
                            if ad.edge[m][m_].snr[s,ny,nx] > edgethr \
                               and ad.cross[(m,m_)].snr[s,ny,nx] > crossthr:
                                idx[s,ny,nx] = p
                                p += 1
                idxt.append(idx)
            edge.append(idxt)
        return (edge, p)
        
    def makeintraindex(self, ad, p, montthr):
        intra = []
        if True: #ad.S>1:
            for m in range(ad.M):
                S,NY,NX = ad.intra[m].shape()
                idx = np.zeros((S,NY,NX), dtype=int) - 1
                for ny in range(NY):
                    for nx in range(NX):
                        for s in range(S):
                            if ad.intra[m].snr[s,ny,nx]  > montthr:
                                idx[s,ny,nx] = p
                                p += 1
                intra.append(idx)
        return (intra, p)
        
class Matrix:
    def __init__(self, ad, idx, dim='x',
                 w_cross=1, w_intra=1, w_edge=1, w_elast=.01):
        # ad must be of type AllDeltas
        # idx must be of type Index
        # dim must be either "x" or "y".
        self.idx = idx
        self.P = idx.P
        self.diag = np.zeros(self.P)
        self.upper = np.zeros(self.P-1)
        self.lower = np.zeros(self.P-1)
        self.b = np.zeros(self.P)
        self.add_E_elast(ad, dim, w_elast)
        self.add_E_cross(ad, dim, w_cross)
        self.add_E_intra(ad, dim, w_intra)
        self.add_E_edge(ad, dim, w_edge)
        self.A = scipy.sparse.diags([self.diag,self.upper,self.lower], \
                                    [0,1,-1], format='csr')

    def add_E_elast(self, ad, dim, alpha):
        for m in range(ad.M):
            S,NY,NX = ad.intra[m].shape()
            for s in range(S):
                for ny in range(NY):
                    for nx in range(NX):
                        idx = self.idx.intra[m][s,ny,nx]
                        if idx>=0:
                            self.diag[idx] += 2*alpha
            # Need to add elastic constraints to edge points also?
            
    def add_E_cross(self, ad, dim, gamma):
        for m in range(ad.M):
            for m_ in range(ad.M):
                mm = (m,m_)
                self.add_one_cross(ad.cross[(m,m_)], self.idx.edge[m][m_],
                                   ad.montpos[m], ad.montpos[m_],
                                   dim, gamma)

    def add_one_cross(self, dcrs, icrs, montpos, montpos2, dim, gamma):
                S,NY,NX = dcrs.shape()
                for s in range(S):
                    for ny in range(NY):
                        for nx in range(NX):
                            idx = icrs[s,ny,nx]
                            if idx>=0:
                                self.diag[idx] += 2*gamma # for X
                                if dim=='x':
                                    dmont = montpos[0] - montpos2[0]
                                    dd = dcrs.dx[s,ny,nx] - dmont/2
                                else:
                                    dmont = montpos[1] - montpos2[1]
                                    dd = dcrs.dy[s,ny,nx] - dmont/2
                                self.b[idx] += 2*gamma * dd

    def add_E_intra(self, ad, dim,  beta):
        for m in range(ad.M):
            self.add_one_intra(ad.intra[m], self.idx.intra[m], dim, beta)
            
    def add_one_intra(self, dintra, iintra, dim, beta):
        S,NY,NX = dintra.shape()
        for s in range(1,S):
            for ny in range(NY):
                for nx in range(NX):
                    idx = iintra[s,ny,nx]
                    idx1 = iintra[s-1,ny,nx]
                    if idx>=0 and idx1>0:
                        if idx1 != idx-1:
                            raise Exception(f'Expected idx to match at' \
                                            + f'{m},{k},{s}')
                        self.diag[idx] += 2*beta
                        self.diag[idx1] += 2*beta
                        self.upper[idx1] += -2*beta
                        self.lower[idx1] += -2*beta
                        if dim=='x':
                            dd = dintra.dx[s,ny,nx] - dintra.dx[s-1,ny,nx]
                        else:
                            dd = dintra.dy[s,ny,nx] - dintra.dy[s-1,ny,nx]
                        self.b[idx] += 2*beta * dd
                        self.b[idx1] -= 2*beta * dd

    def add_E_edge(self, ad, dim, beta):
        for m in range(ad.M):
            for m_ in range(ad.M):
                self.add_one_edge(ad.edge[m][m_], self.idx.edge[m][m_],
                                  dim, beta)
                
    def add_one_edge(self, dedge, iedge, dim, beta):
        S,NY,NX = dedge.shape()
        for s in range(1,S):
            for ny in range(NY):
                for nx in range(NX):
                    idx = iedge[s,ny,nx]
                    idx1 = iedge[s-1,ny,nx]
                    if idx>=0 and idx1>=0:
                        self.diag[idx] += 2*beta
                        self.diag[idx1] += 2*beta
                        self.upper[idx1] += -2*beta
                        self.lower[idx1] += -2*beta
                        if dim=='x':
                            dd = dedge.dx[s,ny,nx] - dedge.dx[s-1,ny,nx]
                        else:
                            dd = dedge.dy[s,ny,nx] - dedge.dy[s-1,ny,nx]
                        self.b[idx] += 2*beta * dd
                        self.b[idx1] -= 2*beta * dd

class Solution:
    def __init__(self, ad, matx, maty):
        dx = scipy.sparse.linalg.spsolve(matx.A, matx.b)
        dy = scipy.sparse.linalg.spsolve(maty.A, maty.b)
        self.intra = self.extractintra(ad.intra, matx.idx.intra, dx, dy)
        self.edge = self.extractedge(ad.edge, matx.idx.edge, dx, dy)
        
    def extractintra(self, adintra, idxintra, dx, dy):
        dlts = []
        M = len(adintra)
        for m in range(M):
            S,NY,NX = adintra[m].shape()
            dlt = Delta((S,NY,NX))
            dlt.ix = adintra[m].ix
            dlt.iy = adintra[m].iy
            dlt.xx = adintra[m].xx
            dlt.yy = adintra[m].yy
            dlt.dx += np.nan
            dlt.dy += np.nan
            for ny in range(NY):
                for nx in range(NX):
                    for s in range(S):
                        idx = idxintra[m][s,ny,nx]
                        if idx>=0:
                            dlt.dx[s,ny,nx] = dx[idx]
                            dlt.dy[s,ny,nx] = dy[idx]
            dlts.append(dlt)
        return dlts
            
    def extractedge(self, adedge, idxedge, dx, dy):
        M = len(adedge)
        dlts = []
        for m in range(M):
            dlts.append(self.extractintra(adedge[m], idxedge[m], dx, dy))
            # Surprisingly, the extractintra code can be used to unpack
            # one set of edges as well.
        return dlts

class UDelta:
    def __init__(self, SHP):
        # SHP is (S, NY, NX), where NY = 2 + (number of vertical tiles)
        # and NX = 2 + (number of horizontal tiles).
        # That's because we include the edge positions. (Red dots in
        # figure on E&R p. 1577.)
        self.xx = np.zeros(SHP) + np.nan
        self.yy = np.zeros(SHP) + np.nan
        self.dx = np.zeros(SHP) + np.nan
        self.dy = np.zeros(SHP) + np.nan

    def infermissing(self):
        # For now, this operates entirely within a tile
        # If anything is missing, we don't trust other data, so:
        self.dy[np.isnan(self.dx)] = np.nan
        self.dx[np.isnan(self.dy)] = np.nan
        self.xx[np.isnan(self.dx)] = np.nan
        self.yy[np.isnan(self.dx)] = np.nan
        self.supported = np.logical_not(np.isnan(self.dx))
        self.infermissingxx()
        self.infermissingyy()
        self.dx = self.infermissing1d(self.dx)
        self.dy = self.infermissing1d(self.dy)

    def infermissingxx(self):
        S,NY,NX = self.xx.shape
        for s in range(S):
            xx = self.xx[s,:,:]
            for ix in range(NX):
                isn = np.isnan(xx[:,ix])
                if np.all(isn):
                    if ix==0:
                        xx[:,ix] = 0
                    elif ix==NX-1:
                        xx[:,ix] = MAXX
                    else:
                        xx[:,ix] = .5*X + (ix-1)
                else:
                    xx[isn, ix] = np.mean(xx[np.logical_not(isn), ix])
            self.xx[s,:,:] = xx

    def infermissingyy(self):
        S,NY,NX = self.yy.shape
        for s in range(S):
            yy = self.yy[s,:,:]
            for iy in range(NY):
                isn = np.isnan(yy[iy,:])
                if np.all(isn):
                    if iy==0:
                        yy[iy,:] = 0
                    elif iy==NY-1:
                        yy[iy,:] = MAXX
                    else:
                        yy[iy,:] = .5*Y + (iy-1)
                else:
                    yy[iy, isn] = np.mean(yy[iy, np.logical_not(isn)])
            self.yy[s,:,:] = yy

    def infermissing1d(self, dd):
        S,NY,NX = self.xx.shape
        for s in range(S):
            d0 = dd[s,:,:]
            d1 = d0.copy() 
            for ny in range(NY):
                for nx in range(NX):
                    if np.isnan(d0[ny,nx]):
                        dist2 = (self.xx[s,:,:] - self.xx[s,ny,nx])**2 \
                                + (self.yy[s,:,:] - self.yy[s,ny,nx])**2
                        dist2[np.isnan(d0)] = np.inf
                        dist2 /= np.min(dist2)
                        wei = np.exp(-dist2)
                        use = np.logical_not(np.isnan(d0))
                        d1[ny,nx] = np.sum(wei[use]*d0[use]) / (np.sum(wei[use])+1e-99)
            dd[s,:,:] = d1
        return dd
    
class UnifiedSolution:
    # Only relevant member is "deltas" which is an M-long list of UDeltas.
    # Use infermissing() method to replace nans with reasonable estimates
    # based on neighbors.
    def __init__(self, soln):
        self.deltas = []
        M = len(soln.intra)
        for m in range(M):
            # start with internal points
            dlt = self.makeintradelta(soln.intra[m])

            # now add edge points
            for m_ in range(M):
                if not soln.edge[m][m_].empty():
                    dlt = self.addedgedelta(dlt, soln.edge[m][m_], m<m_)

            self.deltas.append(dlt)

    def infermissing(self):
        # For now, this operates entirely within a tile
        for dlt in self.deltas:
            dlt.infermissing()
            
    def makeintradelta(self, soln):
        (S,NY,NX) = soln.shape()
        dlt = UDelta((S, NY+2, NX+2))
        dlt.xx[:,1:-1,1:-1] = soln.xx
        dlt.yy[:,1:-1,1:-1] = soln.yy
        dlt.dx[:,1:-1,1:-1] = soln.dx
        dlt.dy[:,1:-1,1:-1] = soln.dy
        return dlt

    def addedgedelta(self, dlt, soln, mlessm):
        (S,NY,NX) = soln.shape()
        if NX==1:
            # Side-by-side overlap
            if mlessm:
                # We are the left tile
                myidx = -1
            else:
                # We are the right tile
                myidx = 0
            dlt.xx[:,1:-1,[myidx]] = soln.xx
            dlt.yy[:,1:-1,[myidx]] = soln.yy
            dlt.dx[:,1:-1,[myidx]] = soln.dx
            dlt.dy[:,1:-1,[myidx]] = soln.dy
        elif NY==1:
            # Top-to-bottom overlap
            if mlessm:
                # We are the top tile
                myidx = -1
            else:
                # We are the bottom tile
                myidx = 0
            dlt.xx[:,[myidx],1:-1] = soln.xx
            dlt.yy[:,[myidx],1:-1] = soln.yy
            dlt.dx[:,[myidx],1:-1] = soln.dx
            dlt.dy[:,[myidx],1:-1] = soln.dy
        return dlt
