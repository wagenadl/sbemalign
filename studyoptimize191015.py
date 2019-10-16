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
monttbl = 'montagealignq5relhp'
slictbl = 'slicealignq5'
touchtbl = 'montagealignattouchq5relhp'
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
       L = db.sel(f'''select count(*) from {slictbl}
           where r={r} and m1={m} and m2={m_} and s=0''')[0][0]
       LL[m,m_] = LL[m_,m] = L
LLL = np.sum(LL, 1) # Total number of overlap points in each montage
PP = KK + LLL # Total number of points (internal + overlap) in each montage
P = np.sum(PP) # Grant total number of points per slice
# Note that we are treating corresponding points in different montages
# separately for now

# Find out locations and shifts between montages
xstar = {}
ystar = {}
deltaxstar = {}
deltaystar = {}
weights = {}
print('Getting intra-slice shifts from db')
for m in range(M):
    for m_ in range(m+1, M):
        mm_ = (m,m_)
        m_m = (m_,m)
        L = LL[m,m_]
        if L==0:
            next
        (s, ix1,iy1, ix2,iy2, dx0, dy0, dx,dy, snr) = db.vsel(f'''select
            s, ix1,iy1, ix2,iy2, dx0, dy0, dx+dxb+dxc,dy+dyb+dyc, snrc
            from {slictbl}
            where r={r} and m1={m} and m2={m_} order by ix1,iy1,s''')
        if len(s) != S*L:
            raise Exception(f'Mismatched overlap point count M{m}:{m_}')
        x1 = ix1*X + X/2 + dx0/2 - dx/2
        x2 = ix2*X + X/2 - dx0/2 + dx/2
        y1 = iy1*Y + Y/2 + dy0/2 - dy/2
        y2 = iy2*Y + Y/2 - dy0/2 + dy/2
        # These positions are in the full tile
        dx = x2 - x1
        dy = y2 - y1
        SIZ = (L, S)
        xstar[mm_] = np.reshape(x1, SIZ)
        xstar[m_m] = np.reshape(x2, SIZ)
        deltaxstar[mm_] = np.reshape(dx/2, SIZ)
        deltaxstar[m_m] = np.reshape(-dx/2, SIZ)
        ystar[mm_] = np.reshape(y1, SIZ)
        ystar[m_m] = np.reshape(y2, SIZ)
        deltaystar[mm_] = np.reshape(dy/2, SIZ)
        deltaystar[m_m] = np.reshape(-dy/2, SIZ)
        wei = (snr > 20).astype('float') # Hard threshold
        weights[mm_] = np.reshape(wei, SIZ)
        weights[m_m] = np.reshape(wei, SIZ)
        
# Now deltastar[m][m'][k,s] is the kth overlap point between montages m and m'
# in slice s

Delta00 = np.zeros((M,M,2)) # 3rd dim is x,y
for m in range(M):
    for m_ in range(M):
        mm_ = (m,m_)
        if LL[m,m_]>0:
            Delta00[m,m_,0] = np.mean(deltaxstar[mm_])
            Delta00[m,m_,1] = np.mean(deltaystar[mm_])
            # This could be weighted by snr
Delta0 = np.zeros((M,2)) # 2nd dimension is x,y.
for a in range(2):
    Delta0[:,a] = np.sum(Delta00[:,:,a]*LL, 1) / np.sum(LL, 1)
# Delta0 is montage shift wrt model space

# Find out locations and shifts within montages - internal points
print('Getting intra-montage shifts from db [internal]')
xxk = {}
yyk = {}
deltaxk = {}
deltayk = {}
weightk = {}
for m in range(M):
    (s, ix,iy, x,y, dx,dy, snr) = db.vsel(f'''select
        s, ix,iy, -dx, -dy, dx+dxb,dy+dyb, snrb
        from {monttbl}
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
    wei = (snr > 20).astype('float') # Hard threshold
    weightk[m] = np.reshape(wei, SIZ)

# Find out locations and shifts within montages - touch points
print('Getting intra-montage shifts from db [touch]')
xxt = {}
yyt = {}
deltaxt = {}
deltayt = {}
weightt = {}
for m in range(M):
    for m_ in range(M):
        L = LL[m,m_]
        mm_ = (m,m_)
        if L>0:
            (s, ix,iy, x,y, dx,dy, snrb) = db.vsel(f'''select
            s, ix,iy, -dx, -dy, dx+dxb,dy+dyb, snrb from {touchtbl}
            where r={r} and m={m} and m2={m_} order by ix,iy,s''')
            if len(s) != S*L:
                raise Exception(f'Mismatched point count for touch tbl M{m}')
            x += X*ix + X/2
            y += Y*iy + Y/2
            SIZ = (L, S)
            xxt[mm_] = np.reshape(x, SIZ)
            yyt[mm_] = np.reshape(y, SIZ)
            deltaxt[mm_] = np.reshape(dx, SIZ)
            deltayt[mm_] = np.reshape(dy, SIZ)
            wei = (snr > 20).astype('float') # Hard threshold
            weightt[mm_] = np.reshape(wei, SIZ)

# xxt,yyt should be more or less equal to xstar,ystar

# Following are reasonable estimates. I think we don't actually need them.
ptsx = []
ptsy = []
Deltax = []
Deltay = []

# First, let's put in the internal points
print('Estimating shifts for internal points')
idx_internal = np.zeros((M,K,S), dtype=int) - 1
if S>1:
    for m in range(M):
        for k in range(KK[m]):
            for s in range(S):
                if weightk[m][k,s] > 0:
                    idx_internal[m,k,s] = len(ptsx)
                    ptsx.append(xxk[m][k])
                    ptsy.append(xxy[m][k])
                    if s==0:
                        Deltax.append(.5*(deltaxk[m][k,0] - deltaxk[m][k,1]))
                        Deltay.append(.5*(deltayk[m][k,0] - deltayk[m][k,1]))
                    elif s==S-1:
                        Deltax.append(.5*(deltaxk[m][k,-1] - deltaxk[m][k,-2]))
                        Deltay.append(.5*(deltayk[m][k,-1] - deltayk[m][k,-2]))
                    else:
                        Deltax.append(deltaxk[m][k,s] - .5*(deltaxk[m][k,s-1]
                                                            +deltaxk[m][k,s+1]))
                        Deltay.append(deltayk[m][k,s] - .5*(deltayk[m][k,s-1]
                                                            +deltayk[m][k,s+1]))

# Now, let's put in touch points
print('Estimating shifts for touch points')
idx_touch = {}
for m in range(M):
    for m_ in range(M):
        L = LL[m,m_]
        idx = np.zeros((M,L,S), dtype=int)  - 1
        for l in range(L):
            for s in range(S):
                if weightt[mm_][l,s] + weights[mm_][l,s] > 0:
                    mm_ = (m,m_)
                    idx[m,l,s] = len(ptsx)
                    ptsx.append(xstar[mm_][l,s])
                    ptsy.append(ystar[mm_][l,s])
                    delx = deltaxstar[mm_] - Delta0[m,0]
                    dely = deltaystar[mm_] - Delta0[m,1]
                    if S > 1:
                        if s==0:
                            Delx = .5*(deltaxt[mm_][k,0] - deltaxt[mm_][k,1])
                            Dely = .5*(deltayt[mm_][k,0] - deltayt[mm_][k,1])
                        elif s==S-1:
                            Delx = .5*(deltaxt[mm_][k,-1] - deltaxt[mm_][k,-2])
                            Dely = .5*(deltayt[mm_][k,-1] - deltayt[mm_][k,-2])
                        else:
                            Delx = deltaxt[mm_][k,s] - .5*(deltaxt[mm_][k,s-1]
                                                           +deltaxt[mm_][k,s+1])
                            Dely = deltayt[mm_][k,s] - .5*(deltayt[mm_][k,s-1]
                                                           +deltayt[mm_][k,s+1])
                        delx = kappa*delx + (1-kappa) * Delx
                        dely = kappa*dely + (1-kappa) * Dely
                    Deltax.append(delx)
                    Deltay.append(dely)

P = len(ptsx)
pts = np.stack(np.array(ptsx), np.array(ptsy))
Del = np.stack(np.array(Deltax), np.array(Deltay))
A = np.sparse.dok_matrix((2*P,2*P))
b = np.zeros(2*P)

'''
We will solve A Del - b = 0, where A is constructed to minimize an error
function:

E = E_elast + E_mont + E_slic

where

E_elast = α² Σ_ksm [dx_ksm]^2

E_mont = β² Σ_k,s>0,m [(dx_ksm - dx_k,s-1,m) - (deltaxk_ksm - deltaxk_k,s-1,m)]^2
    + β'² Σ_t,s>0,m [(dx_tsm - dx_t,s-1,m) - (deltaxt_tsm - deltaxt_t,s-1,m)]^2

E_slic = γ² Σ_t,s,m,m' [dx_tsm - deltaxstar_ksmm']^2

Solve by dE/d del_p = 0. That yields P linear equations in the Del_p's.
A_pq is the coefficient of Del_q in the p-th equation.
b_p is the constant in the p-th equation. (Mind the sign.)

We use our idx_internal and idx_touch maps to convert between k,m / t,m / t,m,m'
and p indices.
'''

# Let's do E_elast
print('Constructing matrix - E_elast')
alpha = .01
for m in range(M):
    for k in range(KK[m]):
        for s in range(S):
            idx = idx_internal[m,k,s]
            if idx>=0:
                A[idx,idx] = alpha # for X
                A[idx+P,idx+P] = alpha # for Y

# Next, E_slic
print('Constructing matrix - E_slice')
gamma = 1
for m in range(M):
    for m_ in range(M):
        mm_ = (m,m_)
        for l in range(LL[m,m_]):
            for s in range(S):
                idx = idx_touch[mm_][l,s]
                if idx>=0:
                    A[idx, idx] += gamma # for X
                    A[idx+P, idx+P] += gamma # for Y
                    b[idx] = 2*gamma*deltaxstar[mm_][l,s]
                    b[idx+P] = 2*gamma*deltaystar[mm_][l,s]

# Finally, E_mont
# First, the internal points
print('Constructing matrix - E_mont,internal')
beta = 1
for m in range(M):
    for k in range(KK[m]):
        for s in range(1,S):
            idx = idx_internal[m,k,s]
            idx1 = idx_internal[m,k,s-1]
            if idx>=0 and idx1>0:
                A[idx, idx] += beta # for X
                A[idx1, idx1] += beta
                A[idx, idx1] = -beta
                A[idx1, idx] = -beta
                A[idx+P, idx+P] += beta # for Y
                A[idx1+P, idx1+P] += beta
                A[idx+P, idx1+P] = -beta
                A[idx1+P, idx+P] = -beta
                b[idx] += beta * (deltaxk[m][k,s] - deltaxk[m][k,s-1]) # for X
                b[idx1] -= beta * (deltaxk[m][k,s] - deltaxk[m][k,s-1])
                b[idx+P] += beta * (deltayk[m][k,s] - deltayk[m][k,s-1]) # for Y
                b[idx1+P] -= beta * (deltayk[m][k,s] - deltayk[m][k,s-1])

# Then the touch points:
print('Constructing matrix - E_mont,touch')
betap = 1
for m in range(M):
    for m_ in range(M_):
        for k in range(LL[m,m_]):
            for s in range(1,S):
                mm_ = (m,m_)
                idx = idx_touch[mm_][k,s]
                idx1 = idx_touch[mm_][k,s-1]
                if idx>=0 and idx1>0:
                    A[idx, idx] += betap # for X
                    A[idx1, idx1] += betap
                    A[idx, idx1] = -betap
                    A[idx1, idx] = -betap
                    A[idx+P, idx+P] += betap # for Y
                    A[idx1+P, idx1+P] += betap
                    A[idx+P, idx1+P] = -betap
                    A[idx1+P, idx+P] = -betap
                    b[idx] += betap * (deltaxt[mm_][k,s]
                                       - deltaxt[mm_][k,s-1]) # for X
                    b[idx1] -= betap * (deltaxt[mm_][k,s]
                                        - deltaxt[mm_][k,s-1])
                    b[idx+P] += betap * (deltayt[mm_][k,s]
                                         - deltayt[mm_][k,s-1]) # for Y
                    b[idx1+P] -= betap * (deltayt[mm_][k,s]
                                          - deltayt[mm_][k,s-1])

# Seriously. Is that it?
print('Solving matrix')
A1 = scipy.sparse.csr_matrix(A)
Delxy = scipy.sparse.linalg.spsolve(A, b)
Delx = Delxy[:P]
Dely = Delxy[P:]

print('Interpreting results - intra-slice')
# Let's see how well we are doing by comparing our results with desiderata
Delxstar = {}
Delystar = {}
for m in range(M):
    for m_ in range(m+1, M):
        mm_ = (m,m_)
        L = LL[m,m_]
        if L>0:
            Delxstar[mm_] = np.zeros(L,S) + np.nan
            Delystar[mm_] = np.zeros(L,S) + np.nan
            for l in range(L):
                for s in range(S):
                    idx = idx_touch[mm_][l,s]
                    if idx>=0:
                        Delxstar[mm_][l,s] = Delx[idx]
                        Delystar[mm_][l,s] = Dely[idx]
                        
