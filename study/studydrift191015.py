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
r = 1
m = 0
(s,dx,dy) = db.vsel(f'''select s,dx+dxb,dy+dyb from montagealignq25co
                    where r={r} and m={m} and s<700 order by s''')
(sdi,dxdi,dydi) = db.vsel(f'''select s,dx+dxb,dy+dyb from montagealignq25codi
                    where r={r} and m={m} and s<700 order by s''')
S = np.max(s) + 1
s0 = S//2
dxf = dx.copy()
dyf = dy.copy()
LAMBDA = 30
fac = np.exp(-1/LAMBDA)
for s in range(s0+1, S):
    dxf[s] = fac * (dxf[s-1] + dx[s] - dx[s-1])
    dyf[s] = fac * (dyf[s-1] + dy[s] - dy[s-1])
for s in range(s0-1,-1,-1):
    dxf[s] = fac * (dxf[s+1] + dx[s] - dx[s+1])
    dyf[s] = fac * (dyf[s+1] + dy[s] - dy[s+1])

qp.figure('/tmp/s2')
qp.pen('b', 1.5)
qp.plot(np.arange(S), dx)
qp.pen('r', 1.5)
qp.plot(np.arange(S), dy)

qp.pen('005', 1.5)
qp.plot(np.arange(S), dxdi)
qp.pen('500', 1.5)
qp.plot(np.arange(S), dydi)

qp.pen('459', 1.5)
qp.plot(np.arange(S), dxf)
qp.pen('954', 1.5)
qp.plot(np.arange(S), dyf)
qp.pen('k', .5)
qp.xaxis(ticks=np.arange(0,701,100), y=-40)
qp.yaxis(ticks=np.arange(-40,61,20))
qp.plot([0,700], [0,0])
qp.shrink()

