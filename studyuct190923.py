#!/usr/bin/python3

import aligndb
import pyqplot as qp
import numpy as np
import em170428.runinfo
import rawimage

db = aligndb.DB()
r,m,s,dx,dy,sx,sy,snr,iter = db.vsel('select * from uctshift')

qp.figure('/tmp/s1')
qp.marker(size=2)
qp.pen('r')
qp.mark(dx[snr<5], dy[snr<5])
qp.pen('339')
qp.mark(dx[np.logical_and(snr>=5, snr<10)],
        dy[np.logical_and(snr>=5, snr<10)])
qp.pen('k')
qp.mark(dx[snr>=10], dy[snr>=10])

qp.pen('k')
qp.xaxis('dx', ticks=np.arange(-5000, 5000, 1000), y=-10000)
qp.yaxis('dy', ticks=np.arange(-8000, 8000, 2000), x=-5000)
qp.shrink()

qp.figure('/tmp/s2')
qp.marker(size=2)
qp.mark(iter + np.random.rand(len(iter))*.9, snr)
qp.xaxis('iter')
qp.yaxis('snr')
qp.shrink()

ri = em170428.runinfo.RunInfo()
z = ri.runSliceToZ(9, 252)
print(z, z * .050, ri.montageCount(9))

img = rawimage.scaledraw(9, 4, 252)

qp.figure('/tmp/s3', 3, 3)
qp.imsc(img)

im1 = rawimage.betaimg(2310,7)
Y,X = im1.shape
qp.figure('/tmp/s4',X/100,Y/100)
qp.imsc(im1)
