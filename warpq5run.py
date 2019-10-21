#!/usr/bin/python3

import aligndb
import time
import sys
import traceback

import rawimage
import warp
import pyqplot as qp
import numpy as np
import factory

db = aligndb.DB()
ri = db.runinfo()
X = Y = 684*5 # Full size of a q5 image
NY = NX = 7 # Number of measurement points

root = '/lsi2/dw/170428/runalignq5'

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
            xc[m4,-1,-1] = meanx - xm[m4]
            meany = (y[m1,-1,-1]+ym[m1] + y[m2,-1,0]+ym[m2]
                     + y[m3,0,-1]+ym[m3] + y[m4,0,0]+ym[m4]) / 4
            yc[m1,-1,-1] = meany-ym[m1]
            yc[m2,-1,0] = meany - ym[m2]
            yc[m3,0,-1] = meany - ym[m3]
            yc[m4,-1,-1] = meany - ym[m4]
    return (xc, yc)

def roughpos(r):
    (xm, ym) = db.vsel(f'select x, y from roughq5pos where r={r}')
    return (xm, ym)
    
def renderq5quad(mdl, img, x, y, dx, dy, xc, yc, x0, y0, xm, ym, rms):
    # Got positions in order tl, tr, bl, br
    xmodel = x + xm
    ymodel = y + ym
    ximage = x - dx
    yimage = y - dy
    xmdlbox = xc + xm
    ymdlbox = yc + ym
    ovr, msk, xl, yt = warp.warpPerspectiveBoxed(img, xmdlbox, ymdlbox,
                                                 xmodel, ymodel, ximage, yimage)
    try:
        warp.copyWithMask(mdl, ovr, msk, xl-x0, yt-y0)
    except:
        print('Failed to copy', rms)
        raise

def renderq5(mdl, r, m, s, x0, y0, xx, yy, dx, dy, xc, yc, xm, ym):
    print(f'getting images R{r} M{m} S{s}')
    img = fullq5img(r, m, s)
    for nx in range(NX-1):
        for ny in range(NY-1):
            renderq5quad(mdl, img,
                         xx[ny:ny+2,nx:nx+2].flatten(),
                         yy[ny:ny+2,nx:nx+2].flatten(),
                         dx[ny:ny+2,nx:nx+2].flatten(),
                         dy[ny:ny+2,nx:nx+2].flatten(),
                         xc[ny:ny+2,nx:nx+2].flatten(),
                         yc[ny:ny+2,nx:nx+2].flatten(),
                         x0, y0,
                         xm, ym,
                         (r,m,s,nx,ny))

def warpq5run(r, ss=None, usedb=False):
    if ss is None:
        ss = range(0, ri.nslices(r))

    print(f'working on R{r}')
    x0, y0, x1, y1 = runextent(r)
    (xxm, yym) = roughpos(r)

    W = int(x1 - x0 + 1)
    H = int(y1 - y0 + 1)

    ssdone = set()
    if usedb:
        for row in db.sel(f'select s from warpq5rundone where r={r}'):
            ssdone.add(int(row[0]))
    print(f'Done R{r}:', ssdone)
    try:
        os.mkdir(f'{root}/R{r}')
    except:
        pass

    for s in ss:
        if s in ssdone:
            print(f'Skipping R{r} S{s}')
        else:
            print(f'Working on R{r} S{s}')
            (xx, yy, dxx, dyy, m_, nx_, ny_) = measuringpoints(r, s)
            (xxc, yyc) = cornerpoints(xx, yy, xxm, yym, r)
            print(f'got position info R{r} S{s}')

            mdl = np.zeros((H,W), dtype=np.uint8)
            for m in range(ri.nmontages(r)):
                renderq5(mdl, r, m, s, x0, y0,
                     xx[m,:,:], yy[m,:,:],
                     dxx[m,:,:], dyy[m,:,:],
                     xxc[m,:,:], yyc[m,:,:],
                     xxm[m], yym[m])
            cv2.imwrite(f'{root}/R{r}/S{s}.jpg', mdl)
            if usedb:
                db.exe(f'insert into warpq5rundone (r,s) values ({r},{s})')

if __name__ == '__main__':

    def droptable():
        db.exe('''drop table warpq5rundone''')

    def maketable():
        db.exe('''create table if not exists warpq5rundone (
        r integer,
        s integer )''')
    
    import factory
    import cv2

    nthreads = 8
    maketable()
    fac = factory.Factory(nthreads)
    for r0 in range(ri.nruns()):
        r = r0+1
        cnt = db.sel(f'select count(*) from warpq5rundone where r={r}')[0][0]
        if cnt<ri.nslices(r):
            fac.request(warpq5run, r, None, True)
                        
