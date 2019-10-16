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
        #    of intraslice deltas, or the full displacement wrt the stack average
        #    in the case of intramontage deltas
    def shape(self):
        return self.ix.shape
    def empty(self):
        return self.ix.size==0
    def __repr__(self):
        S, NY, NX = self.ix.shape
        return f'Delta[ {S} x {NY} x {NX} ]'

def intraslicedeltas(r, m1, m2, tbl):
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
    '''Given a list of Deltas of one montage relative to all the others,
    calculates the average displacement. Returns a (x, y) pair.'''
    sumx = 0
    sumy = 0
    for delta in deltas:
        sumx += np.mean(delta.dx)
        sumy += np.mean(delta.dy)
    dx = sumx / len(deltas)
    dy = sumy / len(deltas)
    return (dx, dy)

def _montagedeltas(r, where, tbl, name):    
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
    (s, ix,iy, x0,y0, dx,dy, snr) = db.vsel(f'''select
        s, ix,iy, dx, dy, dx+dxb,dy+dyb, snrb
        from {tbl}
        where {where} order by s,iy,ix''')
    if len(s) != S*NY*NX:
        raise Exception(f'Mismatched point count for montage tbl M{m}')
    res.xx = np.reshape(X*ix + X/2 - x0, SHP)
    res.yy = np.reshape(Y*iy + Y/2 - y0, SHP)
    res.dx = np.reshape(dx, SHP)
    res.dy = np.reshape(dy, SHP)
    res.snr = np.reshape(snr, SHP)
    return res

def intramontagedeltas(r, m, tbl):
    '''Finds the points where tiles in montage m in run r where compared
    to their neighbors according to the named database table.
    Returns a Delta structure.
    In the structure, (dx,dy) corresponds to the optimal shift of
    the corresponding point (xx,yy). That is, a pixel in (r,m,s) at
    (xx-dx,yy-dy) corresponds to a pixel at (xx,yy) in montage-global
    coordinates.'''
    return _montagedeltas(r, f'r={r} and m={m}', tbl, f'R{r} M{m}')

def montageedgedeltas(r, m, m2, tbl):
    '''Finds the points where tiles in montage m in run r were compared
    to their neighbors (preceding and following slice) according to the
    named database table. Unlike INTRAMONTAGEDELTAS, this function looks
    specifically for edge points, that is, points connected to a neighboring
    montage. (Called “touch points” elsewhere.)
    Returns a Delta structure.
    In the structure, (dx,dy) corresponds to the optimal shift of
    the corresponding point (xx,yy). That is, a pixel in (r,m,s) at
    (xx-dx,yy-dy) corresponds to a pixel at (xx,yy) in montage-global
    coordinates.
    The shape of the structure should match exactly to the first result
    of INTRASLICEDELTAS(r, m, m2, tbl1) with the appropriate table.
    Note that the results here pertain to points in montage m, not in m2.'''
    return _montagedeltas(r, f'r={r} and m={m} and m2={m2}', tbl,
                          f'R{r} M{m}:{m2}')
