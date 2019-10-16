#!/usr/bin/python3

import aligndb
import time
import sys
import traceback

import swiftir
import pyqplot as qp
import numpy as np
import scipy.sparse
import scipy.sparse.linalg

import optimizing

db = aligndb.DB()
ri = db.runinfo()

r = 1
monttbl = 'montagealignq5relhp'
slictbl = 'slicealignq5'
touchtbl = 'montagealignattouchq5relhp'
kappa = .5

M = ri.nmontages(r)
S = ri.nslices(r)

deltastar = {}
for m in range(M):
    for m_ in range(m,M):
        d1, d2 = optimizing.intraslicedeltas(r, m, m_, slictbl)
        deltastar[(m,m_)] = d1
        deltastar[(m_,m)] = d2
deltak = [] # internal
for m in range(M):
    deltak.append(optimizing.intramontagedeltas(r, m, monttbl))

deltat = [] # touch
for m in range(M):
    delm = []
    for m_ in range(M):
        delm.append(optimizing.montageedgedeltas(r, m, m_, touchtbl))
    deltat.append(delm)

montpos = []
for m in range(M):
    dels = []
    for m_ in range(M):
        dels.append(deltastar[(m,m_)])
    montpos.append(optimizing.montageposition(dels))

# Now, let's construct a linear indexing scheme for all our variables
# We'll put the x coordinates first, followed by all the y coords.

# First, let's put in the internal points
print('Indexing internal points')
idx_internal = {}
p = 0
if S>1:
    for m in range(M):
        S1,NY,NX = deltak[m].shape()
        K = KK[m]
        idx = np.zeros((S1,NY,NX), dtype=int) - 1
        for ny in range(NY):
            for nx in range(NX):
                for s in range(S):
                    if deltak[m].snr[s,ny,nx] > 20:
                        idx[s,ny,nx] = p
                        p += 1
        idx_internal[m] = idx
        
# Now, let's put in touch points
print('Indexing for touch points')
idx_touch = []
for m in range(M):
    idxt = []
    for m_ in range(M):
        mm_ = (m,m_)
        S1,NY,NX = deltat[m][m_].shape()
        idx = np.zeros((S1,NY,NX), dtype=int) - 1
        for ny in range(NY):
            for nx in range(NX):
                for s in range(S):
                    if deltat[m][m_].snr[s,ny,nx] > 20:
                        idx[s,ny,nx] = p
                        p += 1
        idxt.append(idx)
    idx_touch.append(idxt)

'''
We will solve A Del - b = 0, where A is constructed to minimize an error
function:

E = E_elast + E_mont + E_slic

where

E_elast = α Σ_ksm [dx_ksm]^2

E_mont = β Σ_k,s>0,m [(dx_ksm - dx_k,s-1,m) - (deltaxk_ksm - deltaxk_k,s-1,m)]^2
    + β' Σ_t,s>0,m [(dx_tsm - dx_t,s-1,m) - (deltaxt_tsm - deltaxt_t,s-1,m)]^2

E_slic = γ Σ_t,s,m,m' [dx_tsm - deltaxstar_ksmm']^2

Solve by dE/d del_p = 0. That yields P linear equations in the Del_p's.
A_pq is the coefficient of Del_q in the p-th equation.
b_p is the constant in the p-th equation. (Mind the sign.)

For instance, dE_slic/d dx_tsmm = 2gamma [dx_tsmm - deltaxstar_ksmm]
so there is a diagonal entry in A: 2gamma, and an entry in b: 2gamma deltaxstar.
Entries need to be added from the various terms in E.

We use our idx_internal and idx_touch maps to convert between k,m / t,m / t,m,m'
and p indices.
'''

P = p
diag = np.zeros(2*P)
upper = np.zeros(2*P-1)
lower = np.zeros(2*P-1)
b = np.zeros(2*P)

# Let's do E_elast
print('Constructing matrix - E_elast')
alpha = .01

for m in range(M):
    S1,NY,NX = deltak[m].shape()
    for s in range(S):
        for ny in range(NY):
            for nx in range(NX):
                idx = idx_internal[m][s,ny,nx]
                if idx>=0:
                    diag[idx] += 2*alpha # for X
                    diag[idx+P] += 2*alpha # for Y

# Next, E_slic
print('Constructing matrix - E_slice')
gamma = 100
for m in range(M):
    for m_ in range(M):
        mm = (m,m_)
        S1,NY,NX = deltastar[mm].shape()
        for s in range(S):
            for ny in range(NY):
                for nx in range(NX):
                    idx = idx_touch[m][m_][s,ny,nx]
                    if idx>=0:
                        diag[idx] += 2*gamma # for X
                        diag[idx+P] += 2*gamma # for Y
                        b[idx] += 2*gamma*(deltastar[mm].dx[s,ny,nx]
                                           - montpos[m][0])
                        b[idx+P] += 2*gamma*(deltastar[mm].dy[s,ny,nx]
                                           - montpos[m][1])

# Finally, E_mont
# First, the internal points
print('Constructing matrix - E_mont,internal')
beta = 1
for m in range(M):
    S1,NY,NX = deltak[m].shape()
    for s in range(1,S):
        for ny in range(NY):
            for nx in range(NX):
                idx = idx_internal[m][s,ny,nx]
                idx1 = idx_internal[m][s-1,ny,nx]
                if idx>=0 and idx1>0:
                    if idx1 != idx-1:
                        raise Exception(f'Expected idx to match at {m},{k},{s}')
                    diag[idx] += 2*beta # for X
                    diag[idx1] += 2*beta
                    upper[idx1] += -2*beta
                    lower[idx1] += -2*beta
                    diag[idx+P] += 2*beta # for Y
                    diag[idx1+P] += 2*beta
                    upper[idx1+P] += -2*beta
                    lower[idx1+P] += -2*beta
                    difx = deltak[m].dx[s,ny,nx] - deltak[m].dx[s-1,ny,nx]
                    dify = deltak[m].dy[s,ny,nx] - deltak[m].dy[s-1,ny,nx]
                    b[idx] += 2*beta * difx # for X
                    b[idx1] -= 2*beta * difx # for X
                    b[idx+P] += 2*beta * dify
                    b[idx1+P] -= 2*beta * dify

# Then the touch points:
print('Constructing matrix - E_mont,touch')
betap = 1
for m in range(M):
    for m_ in range(M):
        S,NY,NX = deltat[m].shape()
        for s in range(1,S):
            for ny in range(NY):
                for nx in range(NX):
                    idx = idx_touch[m][m_][s,ny,nx]
                    idx1 = idx_touch[m][m_][s-1,ny,nx]
                    if idx>=0 and idx1>=0:
                        diag[idx] += 2*betap # for X
                        diag[idx1] += 2*betap
                        upper[idx1] += -2*betap
                        lower[idx1] += -2*betap
                        diag[idx+P] += 2*betap # for Y
                        diag[idx1+P] += 2*betap
                        upper[idx1+P] += -2*betap
                        lower[idx1+P] += -2*betap
                        difx = deltat[m][m_].dx[s,ny,nx] \
                               - deltat[m][m_].dx[s-1,ny,nx]
                        dify = deltat[m][m_].dy[s,ny,nx] \
                               - deltat[m][m_].dy[s-1,ny,nx]
                        b[idx] += 2*betap * difx
                        b[idx1] -= 2*betap * difx
                        b[idx+P] += 2*betap * dify
                        b[idx1+P] -= 2*betap * dify

# Seriously. Is that it?
A = scipy.sparse.diags([diag,upper,lower],[0,1,-1], format='csr')
print('Solving matrix')
Delxy = scipy.sparse.linalg.spsolve(A, b)
Delx = Delxy[:P]
Dely = Delxy[P:]

print('Unpacking results - edge points')
# Let's see how well we are doing by comparing our results with desiderata
Delxstar = {}
Delystar = {}
for m in range(M):
    for m_ in range(M):
        mm_ = (m,m_)
        L = LL[m,m_]
        if L>0:
            Delxstar[mm_] = np.zeros((L,S)) + np.nan
            Delystar[mm_] = np.zeros((L,S)) + np.nan
            for l in range(L):
                for s in range(S):
                    idx = idx_touch[mm_][l,s]
                    if idx>=0:
                        Delxstar[mm_][l,s] = Delx[idx]
                        Delystar[mm_][l,s] = Dely[idx]
                        
print('Unpacking results - internal points')
Delxint = {}
Delyint = {}
for m in range(M):
    K = KK[m]
    Delxint[m] = np.zeros((K,S)) + np.nan
    Delyint[m] = np.zeros((K,S)) + np.nan
    for k in range(K):
        for s in range(S):
            idx = idx_internal[m][k,s]
            if idx>=0:
                Delxint[m][k,s] = Delx[idx]
                Delyint[m][k,s] = Dely[idx]
    
