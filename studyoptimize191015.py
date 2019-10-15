#!/usr/bin/python3

import aligndb
import time
import sys
import traceback

import swiftir
import pyqplot as qp
import numpy as np

import rawimage
import factory

db = aligndb.DB()
ri = db.runinfo()

r = 1
monttbl = 'montagealignq5coa'
slictbl = 'slicealignq5'
touchtbl = 'montagealignattouchq5'
kappa = .5

M = ri.nmontages(r)
S = ri.nslices(r)
X,Y = 684,684 # Size of images at q5

KK = np.zeros(M, dtype=int) # Number of internal points in each montage
for m in range(M):
    KK[m] = db.sel(f'''select count(*) from {monttbl}
    where r={r} and m={m} and s=0''')[0][0]

LL = np.zeros((M,M), dtype=int) # Number of overlap points between pairs
for m in range(M):
    for m_ in range(m+1,M):
        LL[m,m_] = db.sel(f'''select count(*) from {slictbl}
        where r={r} and m1={m} and m2={m_} and s=0''')[0][0]
        LL[m_,m] = LL[m,m_]
LLL = np.sum(LL, 1) # Total number of overlap points in each montage
PP = KK + LLL # Total number of points (internal + overlap) in each montage
P = np.sum(PP) # Grant total number of points per slice
# Note that we are treating corresponding points in different montages
# separately for now

# Find out locations and shifts between montages
xstar = [ [ None for m in range(M) ] for m_ in range(M) ]
ystar = [ [ None for m in range(M) ] for m_ in range(M) ]
deltaxstar = [ [ None for m in range(M) ] for m_ in range(M) ]
deltaystar = [ [ None for m in range(M) ] for m_ in range(M) ]

for m in range(M):
    for m_ in range(m+1,M):
        if LL[m,m_] == 0:
            next
        (s, ix1,iy1, ix2,iy2, dx0, dy0, dx,dy) = db.vsel(f'''select
        s, ix1,iy1, ix2,iy2, dx0, dy0, dx+dxb+dxc,dy+dyb+dyc from {slictbl}
        where r={r} and m1={m} and m2={m_} order by ix1,iy1,s''')
        if len(s) != S*LL[m,m_]:
            raise Exception(f'Mismatched overlap point count M{m}:{m_}')
        x1 = ix1*X + X/2 + dx0/2 - dx/2
        x2 = ix2*X + X/2 - dx0/2 + dx/2
        y1 = iy1*Y + Y/2 + dy0/2 - dy/2
        y2 = iy2*Y + Y/2 - dy0/2 + dy/2
        # These positions are in the full tile
        dx = x2 - x1
        dy = y2 - y1
        SIZ = (LL[m,m_], S)
        xstar[m][m_] = np.reshape(x1, SIZ)
        xstar[m_][m] = np.reshape(x2, SIZ)
        deltaxstar[m][m_] = np.reshape(dx/2, SIZ)
        deltaxstar[m_][m] = np.reshape(-dx/2, SIZ)
        ystar[m][m_] = np.reshape(y1, SIZ)
        ystar[m_][m] = np.reshape(y2, SIZ)
        deltaystar[m][m_] = np.reshape(dy/2, SIZ)
        deltaystar[m_][m] = np.reshape(-dy/2, SIZ)
        
# Now deltastar[m][m'][k,s] is the kth overlap point between montages m and m'
# in slice s

Delta00 = np.zeros((M,M,2)) # 3rd dim is x,y
for m in range(M):
    for m_ in range(M):
        if LL[m,m_] > 0:
            Delta00[m,m_,0] = np.mean(deltaxstar[m][m_])
            Delta00[m,m_,1] = np.mean(deltaystar[m][m_])
            # This could be weighted by snr
Delta0 = np.zeros((M,2)) # 2nd dimension is x,y
for a in range(2):
    Delta0[:,a] = np.sum(Delta00[:,:,a]*LL, 1) / np.sum(LL, 1)

# Find out locations and shifts within montages - internal points
xxk = [ None for m in range(M) ]
yyk = [ None for m in range(M) ]
deltaxk = [ None for m in range(M) ]
deltayk = [ None for m in range(M) ]
for m in range(M):
    (s, ix,iy, x,y, dx,dy) = db.vsel(f'''select
        s, ix,iy, -dx, -dy, dx+dxb,dy+dyb from {monttbl}
        where r={r} and m={m} order by ix,iy,s''')
    if len(s) != S*KK[m]:
        raise Exception(f'Mismatched point count for montage tbl M{m}')
    x += X*ix + X/2
    y += Y*iy + Y/2
    SIZ = (KK[m], S)
    xxk[m] = np.reshape(x, SIZ)
    yyk[m] = np.reshape(y, SIZ)
    deltaxk[m] = np.reshape(dx, SIZ)
    deltayk[m] = np.reshape(dy, SIZ)

# Find out locations and shifts within montages - touch points
xxt = [ [ None for m in range(M) ] for m_ in range(M) ]
yyt = [ [ None for m in range(M) ] for m_ in range(M) ]
deltaxt = [ [ None for m in range(M) ] for m_ in range(M) ]
deltayt = [ [ None for m in range(M) ] for m_ in range(M) ]
for m in range(M):
    for m_ in range(M):
        if LL[m,m_]>0:
            (s, ix,iy, x,y, dx,dy) = db.vsel(f'''select
            s, ix,iy, -dx, -dy, dx+dxb,dy+dyb from {touchtbl}
            where r={r} and m={m} order by ix,iy,s''')
            if len(s) != S*LLL[m]:
                raise Exception(f'Mismatched point count for touch tbl M{m}')
            x += X*ix + X/2
            y += Y*iy + Y/2
            SIZ = (LLL[m], S)
            xxt[m][m_] = np.reshape(x, SIZ)
            yyt[m][m_] = np.reshape(y, SIZ)
            deltaxt[m][m_] = np.reshape(dx, SIZ)
            deltayt[m][m_] = np.reshape(dy, SIZ)
# xxt,yyt should be more or less equal to xstar,ystar

pts = np.zeros((P,S,2)) # Comprehensive list of control points
Deltas = np.zeros((P,S,2)) # Comprehensive list of shifts
P0 = 0
PP0m = np.zeros(M, dtype=int)
PP0t = np.zeros((M,M), dtype=int)

# First, let's put in the internal points
for m in range(M):
    PP0m[m] = P0
    pts[P0:P0+KK[m],:,0] = xxk[m]
    pts[P0:P0+KK[m],:,1] = yyk[m]
    if S>1:
        Deltas[P0:P0+KK[m],0,0] = .5*(deltaxk[m][:,0] - deltaxk[m][:,1])
        Deltas[P0:P0+KK[m],-1,0] = .5*(deltaxk[m][:,-1] - deltaxk[m][:,-2])
        Deltas[P0:P0+KK[m],1:-1,0] = deltaxk[m][:,1:-1] - .5*(deltaxk[m][:,:-2]
                                                              +deltaxk[m][:,2:])
        Deltas[P0:P0+KK[m],0,1] = .5*(deltayk[m][:,0] - deltayk[m][:,1])
        Deltas[P0:P0+KK[m],-1,1] = .5*(deltayk[m][:,-1] - deltayk[m][:,-2])
        Deltas[P0:P0+KK[m],1:-1,1] = deltayk[m][:,1:-1] - .5*(deltayk[m][:,:-2]
                                                              +deltayk[m][:,2:])
    P0 += KK[m]

# Now, let's put in touch points
for m in range(M):
    for m_ in range(M):
        PP0t[m,m_] = P0
        if LL[m,m_] > 0:
            pts[P0:P0+LL[m,m_],:,0] = xstar[m][m_]
            pts[P0:P0+LL[m,m_],:,1] = ystar[m][m_]
            Deltas[P0:P0+LL[m,m_],:,0] = deltaxstar[m][m_] - Delta0[m,0]
            Deltas[P0:P0+LL[m,m_],:,1] = deltaystar[m][m_] - Delta0[m,1]
            if S>1:
                Del = np.zeros((LL[m,m_], S, 2))
                Del[:,0,0] = .5*(deltaxt[m][m_][:,0] - deltaxt[m][m_][:,1])
                Del[:,-1,0] = .5*(deltaxt[m][m_][:,-1] - deltaxt[m][m_][:,-2])
                Del[:,1:-1,0] = deltaxt[m][m_][:,1:-1] -.5*(deltaxt[m][m_][:,:-2]
                                                        +deltaxt[m][m_][:,2:])
                Del[:,0,1] = .5*(deltayt[m][m_][:,0] - deltayt[m][m_][:,1])
                Del[:,-1,1] = .5*(deltayt[m][m_][:,-1] - deltayt[m][m_][:,-2])
                Del[:,1:-1,1] = deltayt[m][m_][:,1:-1] -.5*(deltayt[m][m_][:,:-2]
                                                        +deltayt[m][m_][:,2:])
                Deltas[P0:P0+LL[m,m_],:,:] *= kappa
                Deltas[P0:P0+LL[m,m_],:,:] += (1-kappa) * Del
            P0 += LL[m,m_]

if P0!=P:
    raise Exception('Error in total number of points')
foos = [] # Comprehensive list of functions to minimize
