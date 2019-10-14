#!/usr/bin/python3

# Why is there a jolt b/w R1 M0 S286 and S287?

import aligndb
import time
import sys
import traceback

import swiftir
import numpy as np
import pyqplot as qp
import rawimage

db = aligndb.DB()
ri = db.runinfo()

r = 1
m = 0

for s0 in range(280, 300):
    ss = [s0,s0+1]
    img = []
    for k in range(2):
        img.append(rawimage.q25img(r,m,ss[k]))
    
    dxx = []
    dyy = []
    for k in range(2):
        s = ss[k]
        dx,dy = db.sel(f'''select dx+dxb,dy+dyb from montagealignq25xco
        where r={r} and m={m} and s={s}''')[0]
        #dx,dy=0,0
        dxx.append(dx)
        dyy.append(dy)
        
    Y,X = img[0].shape
    for k in range(2):
        img[k] = swiftir.extractStraightWindow(img[k], (X/2-dxx[k],Y/2-dyy[k]),
                                               (X//4,Y//4))
    
    Y,X = img[0].shape
    qp.figure('/tmp/s1', 6, 3)
    for k in range(2):
        qp.subplot(1,2,k+1)
        qp.imsc(img[k])
        qp.marker('+', 7)
        qp.pen('r',1)
        qp.mark(X/2,Y/2)
        qp.shrink(1,1)
    
    qp.figure('/tmp/s2', 3, 3)
    rgb = np.stack((img[0],img[1],0*img[0]), 2)
    rgb = rgb - np.min(rgb)
    rgb = rgb / np.max(rgb)
    qp.image(rgb)
    time.sleep(1)
