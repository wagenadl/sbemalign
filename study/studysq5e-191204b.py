#!/usr/bin/python3

import aligndb
import swiftir
import rawimage
import pyqplot as qp
import numpy as np

X = Y = 684 # Size of tile
IX = IY = 5

z0 = 3404
r = 20
m1 = 2
m2 = 3
ix1 = 4
iy1 = 3
ix2 = 0
iy2 = 3
s = 79

db = aligndb.DB()
ri = db.runinfo()

img1 = rawimage.partialq5img(r, m1, s, ix1, iy1)
img2 = rawimage.partialq5img(r, m2, s, ix2, iy2)

xm1, ym1 = db.sel(f'''select x,y from solveq5mont
   where r={r} and m={m1}''')[0]
xm2, ym2 = db.sel(f'''select x,y from solveq5mont
   where r={r} and m={m2}''')[0]
xr1, yr1 = db.sel(f'''select x,y from solveq5rigidtile
   where r={r} and m={m1} and s={s}''')[0]
xr2, yr2 = db.sel(f'''select x,y from solveq5rigidtile
   where r={r} and m={m2} and s={s}''')[0]

xx1, yy1, dx1, dy1 = db.vsel(f'''select x, y, dx, dy from solveq5elastic
   where r={r} and m={m1} and s={s} and z0={z0}''')

xx2, yy2, dx2, dy2 = db.vsel(f'''select x, y, dx, dy from solveq5elastic
   where r={r} and m={m2} and s={s} and z0={z0}''')

ii = np.arange(15)
xx0 = xm1 + ((ii+1)%3 - 1)*(X/5) + (ix1+.7)*X
yy0 = ym1 + (ii+1)*Y/17 + iy1*Y + 50

qp.figure('/tmp/s3', 8, 4)
qp.subplot(1,2,1)
qp.mark(xx1, yy1)
qp.shrink(1,1)
qp.subplot(1,2,2)
qp.mark(xx2, yy2)
qp.shrink(1,1)

qp.figure('/tmp/s1', 8, 4)
qp.subplot(1,2,1)
qp.imsc(img1, xx=np.arange(0,X), yy=-np.arange(0,Y))
qp.marker('o', 8, fill='brush')
qp.pen('r', 1)
qp.mark(xx0 - xm1 - xr1 - ix1*X, -(yy0 - ym1 - yr1 - iy1*Y))
qp.shrink(1,1)

qp.subplot(1,2,2)
qp.imsc(img2, xx=np.arange(0,X), yy=-np.arange(0,Y))
qp.marker('o', 8, fill='brush')
qp.pen('r', 1)
qp.mark(xx0 - xm2 - xr2 - ix2*X, -(yy0 - ym2 - yr2 - iy2*Y))
qp.shrink(1,1)

qp.figure('/tmp/s2', 8, 4)
qp.subplot(1,2,1)
qp.imsc(img1, xx=np.arange(0,X), yy=-np.arange(0,Y))
qp.marker('o', 8, fill='brush')
qp.pen('r', 1)
x1 = xx1 - (xm1 + xr1 + ix1*X)
y1 = yy1 - (ym1 + yr1 + iy1*Y)
use = np.logical_and(np.logical_and(x1>0, x1<X), np.logical_and(y1>0, y1<Y))
qp.mark(x1[use], -y1[use])
qp.marker('o', 1, fill='solid')
x1 += dx1
y1 += dy1
qp.mark(x1[use], -y1[use])
qp.shrink(1,1)

qp.subplot(1,2,2)
qp.imsc(img2, xx=np.arange(0,X), yy=-np.arange(0,Y))
qp.marker('o', 8, fill='brush')
qp.pen('g', 1)
x2 = xx1 - (xm2 + xr2 + ix2*X)
y2 = yy1 - (ym2 + yr2 + iy2*Y)
qp.mark(x2[use], -y2[use])
qp.shrink(1,1)



qp.figure('/tmp/s4', 8, 4)
qp.subplot(1,2,2)
qp.imsc(img2, xx=np.arange(0,X), yy=-np.arange(0,Y))
qp.marker('o', 8, fill='brush')
qp.pen('r', 1)
x2 = xx2 - xm2 - xr2 - ix2*X
y2 = yy2 - ym2 - yr2 - iy2*Y
use = np.logical_and(np.logical_and(x2>0, x2<X), np.logical_and(y2>0, y2<Y))
qp.mark(x2[use], -y2[use])
qp.marker('o', 1, fill='solid')
x2 += dx2
y2 += dy2
qp.mark(x2[use], -y2[use])
qp.shrink(1,1)

