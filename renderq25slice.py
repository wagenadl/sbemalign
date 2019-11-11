#!/usr/bin/python3

intbl = 'solveq5slice'

import aligndb
import sys
import os
import numpy as np
import rawimage
import pyqplot as qp

X = Y = 684*5 # Full tile size in q5 space!

db = aligndb.DB()
ri = db.runinfo()

def renderslice(r, s):
    xm, ym = db.vsel(f'''select x,y from {intbl} 
    where r={r} and s={s} order by m''')
    x0 = np.min(xm)
    y0 = np.min(ym)
    xm -= x0
    ym -= y0
    W5 = int(np.max(xm) + X + 1)
    H5 = int(np.max(ym) + Y + 1)
    W25 = W5//5
    H25 = H5//5
    img = np.zeros((H25, W25), dtype=np.uint8)
    M = len(xm)
    C = ri.ncolumns(r)
    R = ri.nrows(r)
    ims = []
    for m in range(M):
        ims.append(rawimage.q25img(r, m, s))
    for m in range(M):
        xb = int(xm[m]/5)
        yb = int(ym[m]/5)
        img[yb:yb+684,xb:xb+684] = ims[m]
        
    for row in range(R-1, -1, -1):
        for col in range(C-1, -1, -1):
            m = row*C + col
            imm = ims[m]
            if col<C-1:
                mr = row*C + col+1
                x1 = (X + xm[mr] - xm[m])/2
            else:
                x1 = X
            if row<R-1:
                mb = (row+1)*C + col
                y1 = (Y + ym[mb] - ym[m])/2
            else:
                y1 = Y
            if col>0:
                ml = row*C + col-1
                x0 = (X - (xm[m] - xm[ml]))/2
            else:
                x0 = 0
            if row>0:
                mt = (row-1)*C + col
                y0 = (Y - (ym[m] - ym[mt]))/2
            else:
                y0 = 0
            x1 = int(x1/5)
            y1 = int(y1/5)
            x0 = int(x0/5)
            y0 = int(y0/5)
            xb = int(xm[m]/5)
            yb = int(ym[m]/5)
            img[yb+y0:yb+y1, xb+x0:xb+x1] = imm[y0:y1,x0:x1]
    return img

img = renderslice(5, 0)
qp.figure('/tmp/s1', img.shape[1]//200, img.shape[0]//200)
qp.imsc(img)
