#!/usr/bin/python3

# This optimizes the relative positions of all the montages in a subvolume
# It does not optimize the positions of individual points.
# This version is relatively primitive: it only uses cross-montage information
# within runs and global "interrun" displacements, not "transrun montage"
# displacements.

crosstbl = 'slicealignq5'
#intratbl = 'relmontalignq5'
#edgetbl  = 'relmontattouchq5'
intertbl = 'interrunq5'
#transtbl = 'transrunmontq5'

outtbl = 'subvolmontposq5'
IX = IY = 5
X = Y = 684
nz = 200

nthreads = 1

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
from submontage import SubMontage

db = aligndb.DB()
ri = db.runinfo()

def droptable():
    db.nofail(f'drop table {outtbl}')

def createtable():
    db.exe(f'''create table if not exists {outtbl} (
    z0 integer,
    nz integer,
    r integer,
    m integer,
    x float,
    y float )''')


class TransDelta:
    def __init__(self, r):
        m,ix,iy,m2,x,y,snr = db.vsel(f'''select
        m,ix,iy,m2,
        x+dx+dxb,y+dy+dyb,snrb
        from {transtbl}
        where r={r}
        order by m,iy,ix''')
        M = ri.nmontages(r)
        if len(m) != M*IX*IY:
            raise Exception(f'''Miscount from {transtbl} at R{r}:
            {len(m)} != {M}*{IX}*{IY} = {M*IX*IY}''')
        SHP  = (M,IY,IX)
        self.m = np.reshape(m, SHP)
        self.ix = np.reshape(ix, SHP)
        self.iy = np.reshape(iy, SHP)
        self.m2 = np.reshape(m2, SHP)
        self.x = np.reshape(x, SHP)
        self.y = np.reshape(y, SHP)
        self.snr = np.reshape(snr, SHP)
    
class DeltaPool:
    def __init__(self, z0, nz):
        self.submonts = SubMontage.allSubMontages(z0, nz)
        self.getrundelta()
        self.getinterdelta()
        self.integrate()
        self.shift()

    def getrundelta(self):
        self.rundelta = {}
        for sm in self.submonts:
            if sm.r not in self.rundelta:
                self.rundelta[sm.r] = optimizing.AllDeltas(sm.r,
                                                           crosstbl=crosstbl,
                                                           s0=sm.s0,
                                                           s1=sm.s1)
        for r in self.rundelta:
            self.rundelta[r].makemontpos()

    def getinterdelta(self):
        self.rr = [r for r in self.rundelta.keys()]
        self.rr.sort()
        self.interdx = {}
        self.interdy = {}
        for k in range(1, len(self.rr)):
            r = self.rr[k]
            dx,dy,snr = db.sel(f'''select
            dx+dxb,dy+dyb,snrb
            from {intertbl}
            where r2={r}''')[0]
            if snr<20:
                raise Exception(f'Unacceptably low SNR in R{r-1}:R{r}: {snr}')
            self.interdx[r] = dx
            self.interdy[r] = dy

    def integrate(self):
        self.montx = {} # Map from r to list of m
        self.monty = {}
        x0 = 0
        y0 = 0
        for r in self.rr:
            rd = self.rundelta[r]
            mx = []
            my = []
            for m in range(len(rd.montpos)):
                mx.append(rd.montpos[m][0]) 
                my.append(rd.montpos[m][1])
                self.montx[r] = np.array(mx) + x0
                self.monty[r] = np.array(my) + y0
            if r<self.rr[-1]:
                x0 -= self.interdx[r+1]
                y0 -= self.interdy[r+1]                

    def shift(self):
        x0 = np.inf
        y0 = np.inf
        for r in self.rr:
            x0 = np.min((x0, np.min(self.montx[r])))
            y0 = np.min((y0, np.min(self.monty[r])))
        for r in self.rr:
            self.montx[r] -= x0
            self.monty[r] -= y0
        x1 = -np.inf
        y1 = -np.inf
        for r in self.rr:
            x1 = np.max((x1, np.max(self.montx[r])))
            y1 = np.max((y1, np.max(self.monty[r])))
        self.x1 = x1 + X
        self.y1 = y1 + Y
        
def optisub(z0, nz):
    print(f'Working on Z{z0}+{nz}')
    dp = DeltaPool(z0, nz)
    with db.db:
        with db.db.cursor() as c:
            c.execute(f'delete from {outtbl} where z0={z0} and nz={nz}')
            for r in dp.rr:
                for m in range(len(dp.montx[r])):
                    c.execute(f'''insert into {outtbl}
                    (z0,nz,
                    r,m,
                    x,y)
                    values
                    ({z0},{nz},
                    {r},{m},
                    {dp.montx[r][m]}, {dp.monty[r][m]})''')

def perhapsoptisub(z0, nz):
    cnt = db.sel(f'select count(*) from {outtbl} where z0={z0} and nz={nz}')
    if cnt[0][0]==0:
        optisub(z0, nz)

createtable()

fac = factory.Factory(nthreads)
R = ri.nruns()
Z = ri.z0(R) + ri.nslices(R)

for z0 in range(0, Z, nz):
    fac.request(perhapsoptisub, z0, nz)
    
fac.shutdown()
