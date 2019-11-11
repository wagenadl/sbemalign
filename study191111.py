#!/usr/bin/python3

import rawimage
import aligndb
import pyqplot as qp
import numpy as np

db = aligndb.DB()
ri = db.runinfo()

r = 5
m1 = 0
m2 = 2
s = 0
ii = 2

x1,y1,x2,y2,dx,dy = db.sel(f'''select x1,y1,x2,y2,dx+dxb+dxc,dy+dyb+dyc
from slicealignq5 
where r={r} and s={s} and m1={m1} and m2={m2} and ii={ii}''')[0]
if m2==m1+1:
    x1 -= 4 * 684
    y1 -= ii * 684
    y2 -= ii * 684
else:
    x1 -= ii * 684
    x2 -= ii * 684
    y1 -= 4 * 684

img1 = rawimage.partialq5img(r,m1,s,4,ii)
img2 = rawimage.partialq5img(r,m2,s,0,ii)



qp.figure('/tmp/s2', 8, 4)
qp.subplot(1,2,1)
qp.imsc(img1, xx=np.arange(684), yy=-np.arange(684))
qp.marker('+', 5)
qp.pen('g', 2)
qp.mark(x1, -y1)
qp.pen('r')
qp.mark(x1 -dx/2, -(y1-dy/2))
qp.shrink()

qp.subplot(1,2,2)
qp.imsc(img2, xx=np.arange(684), yy=-np.arange(684))
qp.marker('+', 5)
qp.pen('g', 2)
qp.mark(x2, -y2)
qp.pen('r')
qp.mark(x2 + dx/2, -(y2+dy/2))
qp.shrink()

