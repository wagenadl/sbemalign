#!/usr/bin/python3

import aligndb
import time
import sys
import traceback

import rawimage
import warp
import pyqplot as qp
import numpy as np

r = 9
s = 0

db = aligndb.DB()
ri = db.runinfo()
X = Y = 684*5 # Full size of a q5 image
NY = NX = 7 # Number of measurement points

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

def measuringpoints(r, s):
    '''Returns MxNYxNX matrices of x, y, dx, dy and also m, nx, ny'''
    (m, nx, ny, x, y, dx, dy) =  db.vsel(f'''select m, nx, ny, x, y, dx, dy
    from optimizeq5 where r={r} and s={s}
    order by m,ny,nx''')
    SHP = ( ri.nmontages(r), NY, NX)
    m = np.reshape(m, SHP)
    nx = np.reshape(nx, SHP)
    ny = np.reshape(ny, SHP)
    x = np.reshape(x, SHP)
    y = np.reshape(y, SHP)
    dx = np.reshape(dx, SHP)
    dy = np.reshape(dy, SHP)
    return (x, y, dx, dy, m, nx, ny)

def cornerpoints(x, y, xm, ym, r):
    '''x, y must be MxNYxNX position of measuring points
    xm, ym must be M-vectors of rough montage positions.
    r must be run number.
    Returns (xc, yc), being MxNYxNX matrices of corner points.
    Those are the same as the measuring points, except that at corners
    and edges shared between multiple montages, the average of all the
    relevant points is taken.'''
    xc = x + 0.0
    yc = y + 0.0
    NC = ri.ncolumns(r)
    NR = ri.nrows(r)
    # Vertical edges between adjacent tiles
    for ir in range(NR):
        for ic in range(NC-1):
            m1 = ir*NC + ic
            m2 = ir*NC + (ic+1)
            meanx = (x[m1,:,-1]+xm[m1] + x[m2,:,0]+xm[m2]) / 2
            xc[m1,:,-1] = meanx - xm[m1]
            xc[m2,:,0] = meanx - xm[m2]
            meany = (y[m1,:,-1]+ym[m1] + y[m2,:,0]+ym[m2]) / 2
            yc[m1,:,-1] = meany - ym[m1]
            yc[m2,:,0] = meany - ym[m2]
    # Horizontal edges between adjacent tiles
    for ir in range(NR-1):
        for ic in range(NC):
            m1 = ir*NC + ic
            m2 = (ir+1)*NC + ic
            meanx = (x[m1,-1,:]+xm[m1] + x[m2,0,:]+xm[m2]) / 2
            xc[m1,-1,:] = meanx-xm[m1]
            xc[m2,0,:] = meanx-xm[m2]
            meany = (y[m1,-1,:]+ym[m1] + y[m2,0,:]+ym[m2]) / 2
            yc[m1,-1,:] = meany-ym[m1]
            yc[m2,0,:] = meany-ym[m2]
    # Corner pieces
    for ir in range(NR-1):
        for ic in range(NC-1):
            m1 = ir*NC + ic
            m2 = ir*NC + (ic+1)
            m3 = (ir+1)*NC + ic
            m4 = (ir+1)*NC + (ic+1)
            meanx = (x[m1,-1,-1]+xm[m1] + x[m2,-1,0]+xm[m2]
                     + x[m3,0,-1]+xm[m3] + x[m4,0,0]+xm[m4]) / 4
            xc[m1,-1,-1] = meanx-xm[m1]
            xc[m2,-1,0] = meanx - xm[m2]
            xc[m3,0,-1] = meanx - xm[m3]
            xc[m4,0,0] = meanx - xm[m4]
            meany = (y[m1,-1,-1]+ym[m1] + y[m2,-1,0]+ym[m2]
                     + y[m3,0,-1]+ym[m3] + y[m4,0,0]+ym[m4]) / 4
            yc[m1,-1,-1] = meany-ym[m1]
            yc[m2,-1,0] = meany - ym[m2]
            yc[m3,0,-1] = meany - ym[m3]
            yc[m4,0,0] = meany - ym[m4]
    return (xc, yc)

def roughpos(r):
    (xm, ym) = db.vsel(f'select x, y from roughq5pos where r={r}')
    return (xm, ym)
    
def renderq5quad(mdl, img, x, y, dx, dy, xc, yc, x0, y0, xm, ym):
    print(f'renderq5quad r{r} m{m} s{s}')
    # Got positions in order tl, tr, bl, br
    xmodel = x + xm
    ymodel = y + ym
    ximage = x - dx
    yimage = y - dy
    print(xmodel)
    print(ymodel)
    print(ximage)
    print(yimage)
    print(f'x0{x0} y0{y0} xm{xm} ym{ym}')
    print('x = ', x)
    print('y = ', y)
    print('dx = ', dx)
    print('dy = ', dy)
    xmdlbox = xc + xm
    ymdlbox = yc + ym
    ovr, msk, xl, yt = warp.warpPerspectiveBoxed(img, xmdlbox, ymdlbox,
                                                 xmodel, ymodel, ximage, yimage)
    print('ovr.shape', ovr.shape)
    qp.subplot(1,2,1)
    qp.imsc(ovr)
    qp.subplot(1,2,2)
    qp.imsc(msk)
    print('msk.shape', msk.shape)
    print('mdl.shape', mdl.shape)
    warp.copyWithMask(mdl, ovr, msk, xl-x0, yt-y0)

def renderq5(mdl, r, m, s, x0, y0, xx, yy, dx, dy, xc, yc, xm, ym):
    print(f'renderq5 R{r} M{m} S{s}')
    img = fullq5img(r, m, s)
    print(f'got image R{r} M{m} S{s}')
    if m==0:
        nxx = [5]
        nyy = [5]
    elif m==1:
        nxx = [0]
        nyy = [5]
    elif m==2:
        nxx = [5]
        nyy = [0]
    elif m==3:
        nxx = [0]
        nyy = [0]
    elif m==4:
        nxx = [2]
        nyy = [2]
    else:
        nxx = []
        nyy = []
    for nx in nxx:
        for ny in nyy:
            qp.figure(f'/tmp/s{nx}{ny}', 6, 3)
            renderq5quad(mdl, img,
                         xx[ny:ny+2,nx:nx+2].flatten(),
                         yy[ny:ny+2,nx:nx+2].flatten(),
                         dx[ny:ny+2,nx:nx+2].flatten(),
                         dy[ny:ny+2,nx:nx+2].flatten(),
                         xc[ny:ny+2,nx:nx+2].flatten(),
                         yc[ny:ny+2,nx:nx+2].flatten(),
                         x0, y0,
                         xm, ym)

t0 = time.time()
print(f'rendering R{r} S{s}')
x0, y0, x1, y1 = runextent(r)
(xx, yy, dxx, dyy, m_, nx_, ny_) = measuringpoints(r, s)
(xxm, yym) = roughpos(r)
(xxc, yyc) = cornerpoints(xx, yy, xxm, yym, r)

W = int(x1 - x0 + 1)
H = int(y1 - y0 + 1)
print(f'got position info R{r} S{s}')

mdl = np.zeros((H,W), dtype=np.uint8)
for m in range(ri.nmontages(r)):
    renderq5(mdl, r, m, s, x0, y0,
             xx[m,:,:], yy[m,:,:],
             dxx[m,:,:], dyy[m,:,:],
             xxc[m,:,:], yyc[m,:,:],
             xxm[m], yym[m])
t1 = time.time()
print(t1-t0)
