#!/usr/bin/python3

import aligndb
import swiftir
import rawimage
import pyqplot as qp
import numpy as np

X = Y = 684 # Size of tile
IX = IY = 5

z0 = 3504
r = 20
m1 = 0
m2 = 2
ix1 = 2
iy1 = 4
ix2 = 2
iy2 = 0
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
   where r={r} and m={m1} and s={s}''')

xx2, yy2, dx2, dy2 = db.vsel(f'''select x, y, dx, dy from solveq5elastic
   where r={r} and m={m1} and s={s}''')

ii = np.arange(15)
xx0 = xm1 + (ii+1)*X/17 + ix1*X
yy0 = ym1 + ((ii+1)%3 - 1)*(Y/5) + (iy1+.78)*Y

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
