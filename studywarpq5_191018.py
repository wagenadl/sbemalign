#!/usr/bin/python3

import aligndb
import time
import sys
import traceback

import rawimage
import warp
import pyqplot as qp
import numpy as np

r = 1
s = 300

db = aligndb.DB()
ri = db.runinfo()
X = Y = 684*5 # Full size of a q5 image

def runextent(r):
    res = db.sel(f'''select
        min(p.x + dx) as x0, min(p.y + dy) as y0,
        {X}+max(p.x + dx) as x1, {Y}+max(p.y + dy) as y1
        from roughq5pos as p
        inner join optimizeq5 as o on p.r=o.r and p.m=o.m where p.r={r}''')
    x0, y0, x1, y1 = res[0]
    x0 = int(x0)
    y0 = int(y0)
    x1 = int(x1+1)
    y1 = int(y1+1)
    return x0,y0,x1,y1

def fullq5img(r, m, s):
    img = np.zeros((Y,X), dtype=np.uint8)
    for iy in range(5):
        for ix in range(5):
            subimg = rawimage.partialq5img(r, m, s, ix, iy)
            H,W = subimg.shape
            img[H*iy:H*(iy+1), W*ix:W*(ix+1)] = subimg
    return img

def roundedquad(xx, yy):
    # xx,yy must be presented in order tl, tr, bl, br
    # results are presented in order tl, tr, br, bl
    return (np.array([int(xx[0]), int(xx[1]+1), int(xx[3]+1), int(xx[2])]),
            np.array([int(yy[0]), int(yy[1]), int(yy[3]+1), int(yy[2])+1]))

def renderq5quad(mdl, img, r, m, s, nx, ny, x0, y0, xm, ym):
    (nx1, ny1, x, y, dx, dy) = db.vsel(f'''select nx, ny, x, y, dx, dy
    from optimizeq5
    where r={r} and m={m} and s={s}
    and (nx={nx} or nx={nx+1})
    and (ny={ny} or ny={ny+1}) order by ny,nx''')
    # Got positions in order tl, tr, bl, br
    xmodel = x + xm
    ymodel = y + ym
    ximage = x - dx
    yimage = y - dy
    mdl2img = warp.getPerspective(xmodel, ymodel, ximage, yimage)
    ovr, xl, yt = warp.warpPerspective(img, mdl2img) # This is deeply suboptimal
    # .. because we only need a small part of it. But let me first check this.
    xmr, ymr = roundedquad(xmodel, ymodel)
    msk = warp.createClipMask(xmr, ymr, x0, y0, W, H)
    warp.copyWithMask(mdl, ovr, msk, xl-x0, yt-y0)

def renderq5(mdl, r, m, s, x0, y0):
    img = fullq5img(r, m, s)
    (xm, ym) = db.sel(f'select x, y from roughq5pos where r={r} and m={m}')[0]
    for nx in range(6):
        for ny in range(6):
            renderq5quad(mdl, img, r, m, s, nx, ny, x0, y0, xm, ym)
           
x0, y0, x1, y1 = runextent(r)
W = int(x1 - x0 + 1)
H = int(y1 - y0 + 1)

mdl = np.zeros((H,W), dtype=np.uint8)
for m in range(ri.nmontages(r)):
    renderq5(mdl, r, m, s, x0, y0)
