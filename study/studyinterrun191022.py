#!/usr/bin/python3

import aligndb
import time
import sys
import traceback

import rawimage
import warpq5run
import warp
import pyqplot as qp
import numpy as np
import swiftir

db = aligndb.DB()
ri = db.runinfo()

r1 = 21
r2 = 22

img1 = warpq5run.warpedq5img(r1, ri.nslices(r1)-1)
img2 = warpq5run.warpedq5img(r2, 0)

Y1,X1 = img1.shape
Y2,X2 = img2.shape

win1 = swiftir.extractStraightWindow(img1, (X1/2,Y1/2), (2048,4096))
win2 = swiftir.extractStraightWindow(img2, (X2/2,Y2/2), (2048,4096))
apo1 = swiftir.apodize(win1)
apo2 = swiftir.apodize(win2)

qp.figure('/tmp/s1', 4, 8)
qp.imsc(apo1)

qp.figure('/tmp/s2', 4, 8)
qp.imsc(apo2)

(dx, dy, sx, sy, snr) = swiftir.swim(apo1, apo2)

win1 = swiftir.extractStraightWindow(img1, (X1/2-dx/2,Y1/2-dy/2), (2048,4096))
win2 = swiftir.extractStraightWindow(img2, (X2/2+dx/2,Y2/2+dy/2), (2048,4096))
apo1 = swiftir.apodize(win1)
apo2 = swiftir.apodize(win2)

qp.figure('/tmp/s1a', 4, 8)
qp.imsc(apo1)

qp.figure('/tmp/s2a', 4, 8)
qp.imsc(apo2)

(dxb, dyb, sxb, syb, snrb) = swiftir.swim(apo1, apo2)

