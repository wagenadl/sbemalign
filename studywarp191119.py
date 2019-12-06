#!/usr/bin/python3

import rawimage
import renderq5utils
import aligndb
import warp
import numpy as np
import cv2

monttbl = 'solveq5mont'
rigidtbl = 'solveq5rigidtile'
elastbl = 'solveq5elastic'
globtbl = 'q5global'
bboxtbl = 'q5bbox'

X = Y = 684*5 # Full tile size in q5 space!

db = aligndb.DB()
ri = db.runinfo()

r = 5
s = 10

tiles = []
M = ri.nmontages(r)
for m in range(M):
    print(f'Loading montage {m} of {M}')
    tiles.append(rawimage.fullq5img(r, m, s))
    
xl, xr, yt, yb = [], [], [], []
for m in range(M):
    xl1, yt1, xr1, yb1 = renderq5utils.renderlimits(r, m, s)
    xl.append(xl1)
    xr.append(xr1)
    yt.append(yt1)
    yb.append(yb1)    
xl0 = np.min(xl)
xr0 = np.max(xr)
yt0 = np.min(yt)
yb0 = np.max(yb)
W = xr0 - xl0
H = yb0 - yt0

img = np.zeros((H,W), dtype=np.uint8)

for m in range(M):
    xx, yy = renderq5utils.rendergrid(r, m, s)
    IX = len(xx) - 1
    IY = len(yy) - 1
    x0, y0 = renderq5utils.rigidtileposition(r, m, s)
    for ix in range(IX):
        for iy in range(IY):
            print(f'Rendering M{m} IX{ix} IY{iy}')
            xl1 = xx[ix]
            xr1 = xx[ix+1]
            yt1 = yy[iy]
            yb1 = yy[iy+1]
            xywh = [xl1-xl0, yt1-yt0, xr1-xl1, yb1-yt1]
            xmdl = np.array([xx[ix], xx[ix], xx[ix+1], xx[ix+1]])
            ymdl = np.array([yy[iy], yy[iy+1], yy[iy], yy[iy+1]])
            dx, dy = renderq5utils.interpolatedshifts(r, m, s, xmdl, ymdl)
            xtile = xmdl - x0 - dx
            ytile = ymdl - y0 - dy
            warp.warpToRectangle(img, xywh, tiles[m],
                                 xmdl-xl0, ymdl-yt0, xtile, ytile)
            
cv2.imwrite('/tmp/rms1.jpg', img)
