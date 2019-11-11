#!/usr/bin/python3

# This version integrates and high-pass filters the RELMONTALIGNQ5 data
# The **b data is copied verbatim; only the dx,dy are integrated

import aligndb
import time
import sys
import traceback

import numpy as np

db = aligndb.DB()

def droptable():
    db.exe('''drop table montagealignq5relhp''')

def maketable():
    db.exe('''create table if not exists montagealignq5relhp (
    r integer,
    m integer,
    s integer,
    ix integer,
    iy integer,

    dx float,
    dy float,
    sx float,
    sy float,
    snr float,

    dxb float,
    dyb float,
    sxb float,
    syb float,
    snrb float
    )''')

ri = db.runinfo()

def getbaseshifts(r, m, ix, iy):
    s,dx,dy = db.vsel(f'''select s,dx,dy from relmontalignq5
            where r={r} and m={m} and ix={ix} and iy={iy} order by s''')
    S = len(s)
    if S != ri.nslices(r):
        raise Exception(f'Bad number of slices for R{r} M{m} {ix},{iy}')
    s0 = S//2

    dxf = np.zeros(S)
    dyf = np.zeros(S)
    LAMBDA = 25
    fac = np.exp(-1/LAMBDA)
    for s in range(s0+1, S):
        dxf[s] = fac * dxf[s-1] + dx[s]
        dyf[s] = fac * dyf[s-1] + dy[s]
    for s in range(s0-1,-1,-1):
        dxf[s] = fac * dxf[s+1] - dx[s+1]
        dyf[s] = fac * dyf[s+1] - dy[s+1]
    return (dxf, dyf)
    
def alignmanysubtiles(r, m, ix, iy):
    cnt = db.sel(f'''select count(1) from montagealignq5relhp
    where r={r} and m={m} and ix={ix} and iy={iy}''')[0][0]
    if cnt==ri.nslices(r):
        return
    S = ri.nslices(r)
    print(f'Working on R{r} M{m} {ix},{iy} [S={S}]')
    if cnt>0:
        db.exe(f'''delete from montagealignq5relhp 
        where r={r} and m={m} and ix={ix} and iy={iy}''')
    dxf, dyf = getbaseshifts(r, m, ix, iy)
    (sx,sy,snr, dxb,dyb,sxb,syb,snrb) = db.vsel(f'''select
        sx,sy,snr, dxb,dyb,dxb,dyb,snrb
        from  relmontalignq5
        where r={r} and m={m} and ix={ix} and iy={iy} order by s''')
    if len(snr) != S:
        raise Exception(f'Bad number of slices for R{r} M{m} {ix},{iy}')

    with db.db:
        with db.db.cursor() as c:
            for s in range(ri.nslices(r)):
                c.execute(f'''insert into montagealignq5relhp 
                (r,m,s,ix,iy,
                dx,dy,sx,sy,snr,
                dxb,dyb,sxb,syb,snrb)
                values
                ({r},{m},{s},{ix},{iy},
                {dxf[s]},{dyf[s]},{sx[s]},{sy[s]},{snr[s]}, 
                {dxb[s]},{dyb[s]},{sxb[s]},{syb[s]},{snrb[s]})''')

def queuealignmontage(r, m):
    cnt = db.sel(f'''select count(1) from montagealignq5relhp
    where r={r} and m={m}''')[0][0]
    N = ri.nslices(r) * 25
    if cnt==N:
        return
    print(f'Considering run {r} montage {m}: {cnt}<{N}')
    for ix in range(5):
        for iy in range(5):
            alignmanysubtiles(r, m, ix, iy)

maketable()
    
for r0 in range(ri.nruns()):
    r  = r0 + 1
    cnt = db.sel(f'''select count(1) from montagealignq5relhp
    where r={r}''')[0][0]
    N = ri.nmontages(r) * ri.nslices(r) * 25
    if cnt < N:
        print(f'Considering run {r}: {cnt}<{N}')
        for m in range(ri.nmontages(r)):
            queuealignmontage(r, m)

