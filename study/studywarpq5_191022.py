#!/usr/bin/python3

import aligndb
import time
import sys
import traceback

import rawimage
import warp
import pyqplot as qp
import numpy as np
import factory

db = aligndb.DB()
ri = db.runinfo()

r = 27
s = 34
m1 = 9
m2 = 11

(nx1, x1,y1,dx1,dy1, sup1) = db.vsel(f'''select nx, x,y,dx,dy, supported
                          from optimizeq5
                          where r={r} and s={s} and m={m1} and ny=6
                          order by nx''')

(nx2, x2,y2,dx2,dy2, sup2) = db.vsel(f'''select nx, x,y,dx,dy, supported
                          from optimizeq5
                          where r={r} and s={s} and m={m2} and ny=0
                          order by nx''')

(xm1,ym1) = db.vsel(f'''select x,y from roughq5pos where r={r} and m={m1}''')
(xm2,ym2) = db.vsel(f'''select x,y from roughq5pos where r={r} and m={m2}''')

print(y1+ym1)
print(y2+ym2)

(ix1, dx,dy, snr) = db.vsel(f'''select ix1, dx+dxb,dy+dyb,snr
                          from slicealignq5
                          where r={r} and s={s} and m1={m1} and iy1=4
                          order by ix1''')
xc1 = 684*ix1 + 684/2
yc1 = 342 + 342/2 
