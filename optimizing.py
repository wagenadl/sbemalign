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
    '''Given a list of cross Deltas of one montage relative to all the others,
    calculates the average displacement. Returns a (x, y) pair.'''
    sumx = 0
    sumy = 0
    n = 0
    for delta in deltas:
        if delta.dx.size>0:
            sumx += np.mean(delta.dx)
            sumy += np.mean(delta.dy)
            n += 1
    dx = sumx / n
    dy = sumy / n
    return (dx, dy)

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
        raise Exception(f'Mismatched point count for montage tbl M{m}')
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
        self.montpos = []
        for m in range(self.M):
            dels = []
            for m_ in range(self.M):
                dels.append(self.cross[(m,m_)])
            self.montpos.append(montageposition(dels))

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
                               or ad.cross[(m,m_)].snr[s,ny,nx] > crossthr:
                                idx[s,ny,nx] = p
                                p += 1
                idxt.append(idx)
            edge.append(idxt)
        return (edge, p)
        
    def makeintraindex(self, ad, p, montthr):
        intra = []
        if ad.S>1:
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
    def __init__(self, ad, idx, w_cross=1, w_intra=1, w_edge=1, w_elast=.01):
        # ad must be of type AllDeltas
        # idx must be of type Index
        self.idx = idx
        self.P = idx.P
        self.diag = np.zeros(2*self.P)
        self.upper = np.zeros(2*self.P-1)
        self.lower = np.zeros(2*self.P-1)
        self.b = np.zeros(2*self.P)
        self.add_E_elast(ad, w_elast)
        self.add_E_cross(ad, w_cross)
        self.add_E_intra(ad, w_intra)
        self.add_E_edge(ad, w_edge)
        self.A = scipy.sparse.diags([self.diag,self.upper,self.lower], \
                                    [0,1,-1], format='csr')

    def add_E_elast(self, ad, alpha):
        for m in range(ad.M):
            S,NY,NX = ad.intra[m].shape()
            for s in range(S):
                for ny in range(NY):
                    for nx in range(NX):
                        idx = self.idx.intra[m][s,ny,nx]
                        if idx>=0:
                            self.diag[idx] += 2*alpha # for X
                            self.diag[idx+self.P] += 2*alpha # for Y
            # Need to add elastic constraints to edge points also?
            
    def add_E_cross(self, ad, gamma):
        for m in range(ad.M):
            for m_ in range(ad.M):
                mm = (m,m_)
                self.add_one_cross(ad.cross[(m,m_)], self.idx.edge[m][m_],
                                   ad.montpos[m], gamma)

    def add_one_cross(self, dcrs, icrs, montpos, gamma):
                S,NY,NX = dcrs.shape()
                for s in range(S):
                    for ny in range(NY):
                        for nx in range(NX):
                            idx = icrs[s,ny,nx]
                            if idx>=0:
                                self.diag[idx] += 2*gamma # for X
                                self.diag[idx+self.P] += 2*gamma # for Y
                                self.b[idx] += 2*gamma \
                                        * (dcrs.dx[s,ny,nx] - montpos[0])
                                self.b[idx+self.P] += 2*gamma \
                                        * (dcrs.dy[s,ny,nx] -montpos[1])

    def add_E_intra(self, ad, beta):
        for m in range(ad.M):
            self.add_one_intra(ad.intra[m], self.idx.intra[m], beta)
            
    def add_one_intra(self, dintra, iintra, beta):
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
                        self.diag[idx] += 2*beta # for X
                        self.diag[idx1] += 2*beta
                        self.upper[idx1] += -2*beta
                        self.lower[idx1] += -2*beta
                        self.diag[idx+self.P] += 2*beta # for Y
                        self.diag[idx1+self.P] += 2*beta
                        self.upper[idx1+self.P] += -2*beta
                        self.lower[idx1+self.P] += -2*beta
                        difx = dintra.dx[s,ny,nx] \
                               - dintra.dx[s-1,ny,nx]
                        dify = dintra.dy[s,ny,nx] \
                               - dintra.dy[s-1,ny,nx]
                        self.b[idx] += 2*beta * difx
                        self.b[idx1] -= 2*beta * difx
                        self.b[idx+self.P] += 2*beta * dify
                        self.b[idx1+self.P] -= 2*beta * dify

    def add_E_edge(self, ad, beta):
        for m in range(ad.M):
            for m_ in range(ad.M):
                self.add_one_edge(ad.edge[m][m_], self.idx.edge[m][m_], beta)
                
    def add_one_edge(self, dedge, iedge, beta):
        S,NY,NX = dedge.shape()
        for s in range(1,S):
            for ny in range(NY):
                for nx in range(NX):
                    idx = iedge[s,ny,nx]
                    idx1 = iedge[s-1,ny,nx]
                    if idx>=0 and idx1>=0:
                        self.diag[idx] += 2*beta # for X
                        self.diag[idx1] += 2*beta
                        self.upper[idx1] += -2*beta
                        self.lower[idx1] += -2*beta
                        self.diag[idx+self.P] += 2*beta # for Y
                        self.diag[idx1+self.P] += 2*beta
                        self.upper[idx1+self.P] += -2*beta
                        self.lower[idx1+self.P] += -2*beta
                        difx = dedge.dx[s,ny,nx] - dedge.dx[s-1,ny,nx]
                        dify = dedge.dy[s,ny,nx] - dedge.dy[s-1,ny,nx]
                        self.b[idx] += 2*beta * difx
                        self.b[idx1] -= 2*beta * difx
                        self.b[idx+self.P] += 2*beta * dify
                        self.b[idx1+self.P] -= 2*beta * dify

class Solution:
    def __init__(self, ad, mat):
        dxy = scipy.sparse.linalg.spsolve(mat.A, mat.b)
        P = len(mat.b)//2
        dx = dxy[:P]
        dy = dxy[P:]
        self.intra = self.extractintra(ad.intra, mat.idx.intra, dx, dy)
        self.edge = self.extractedge(ad.edge, mat.idx.edge, dx, dy)
        
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
        self.xx = np.zeros(SHP) + np.nan
        self.yy = np.zeros(SHP) + np.nan
        self.dx = np.zeros(SHP) + np.nan
        self.dy = np.zeros(SHP) + np.nan
    
class UnifiedSolution:
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
