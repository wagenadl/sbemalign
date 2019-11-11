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

import optimizingrel as optimizing
import factory

TEST = False

nthreads = 12

crosstbl = 'slicealignq5'
intratbl='relmontalignq5'
edgetbl='relmontattouchq5'
outtbl = 'optimizeq5rel'
roughtbl = 'roughq5posrel'

db = aligndb.DB()
ri = db.runinfo()

if TEST:
    nthreads = 1

def droptable():
    db.nofail(f'drop table {outtbl}')
    db.nofail(f'drop table {roughtbl}')

def maketable():
    db.exe(f'''create table if not exists {outtbl} (
    r integer,
    m integer,
    s integer,
    nx integer,
    ny integer,

    x float,
    y float,
    dx float,
    dy float,
    supported boolean
    )''')

    db.exe(f'''create table if not exists {roughtbl} (
    r integer,
    m integer,
    x float,
    y float
    )''')
    

def insertintodb(r, m, delta, dbcon):
    [S, NY, NX] = delta.xx.shape
    for s in range(S):
        if TEST:
            print(f'Inserting R{r} M{m} S{s}')
        for ny in range(NY):
            for nx in range(NX):
                dbcon.execute(f'''insert into {outtbl}
                (r,m,s,nx,ny,
                x,y,dx,dy,supported)
                values
                ({r},{m},{s},{nx},{ny},
                {delta.xx[s,ny,nx]},
                {delta.yy[s,ny,nx]},
                {delta.dx[s,ny,nx]},
                {delta.dy[s,ny,nx]},
                {delta.supported[s,ny,nx]})
                ''')

def insertintorough(r, m, pos, dbcon):
    print(f'Inserting rough R{r} M{m}')
    dbcon.execute(f'''insert into {roughtbl}
    (r,m,x,y) values ({r},{m},{pos[0]},{pos[1]})''')

def optimizerun(r):
    print(f'Working on R{r}')
    print(f'Gathering deltas R{r}')
    deltas = optimizing.AllDeltas(r,
                                  crosstbl=crosstbl,
                                  intratbl=intratbl,
                                  edgetbl=edgetbl)
    print(f'Calculating overall montage positions R{r}')
    deltas.makemontpos()
    
    print(f'Storing  overall montage positions R{r}')
    with db.db:
        with db.db.cursor() as c:
            c.execute(f'delete from {roughtbl} where r={r}')
            for m in range(len(deltas.montpos)):
                insertintorough(r, m, deltas.montpos[m], c)
    
    print(f'Creating index R{r}')
    idx = optimizing.Index(deltas)
    print(f'Creating matrices R{r}')
    matx = optimizing.Matrix(deltas, idx, 'x')
    maty = optimizing.Matrix(deltas, idx, 'y')
    print(f'Solving matrices R{r}')
    soln = optimizing.Solution(deltas, matx, maty)
    print(f'Collecting results R{r}')
    usol = optimizing.UnifiedSolution(soln)
    usol.infermissing()

    print(f'Inserting into database R{r}')
    with db.db:
        with db.db.cursor() as c:
            c.execute(f'delete from {outtbl} where r={r}')
            for m in range(len(usol.deltas)):
                insertintodb(r, m, usol.deltas[m], c)

maketable()
fac = factory.Factory(nthreads)
rr,cnts = db.vsel(f'select r,count(*) from {edgetbl} group by r order by r')
rrr = set()
for r in rr:
    rrr.add(r)
for r0 in range(ri.nruns()):
    r = r0 + 1
    if r in rrr:
        cnt = db.sel(f'select count(*) from {outtbl} where r={r}')[0][0]
        if cnt < 7*7*ri.nmontages(r)*ri.nslices(r):
            fac.request(optimizerun, r)

print('waiting for completion')
fac.shutdown()
 
