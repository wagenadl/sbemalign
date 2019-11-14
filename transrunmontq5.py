#!/usr/bin/python3

import aligndb
import time
import sys
import traceback

import rawimage
import pyqplot as qp
import numpy as np
import swiftir
import factory

db = aligndb.DB()
ri = db.runinfo()

roughtbl = 'solveq5slice'
intbl = 'interrunq25'
inq = 25
outq = 5
outtbl= 'transrunmontq5'

X = Y = 684 # Size of tile
IX = IY = 5
PICTURES = False
nthreads = 12

def maketable():
    db.exe(f'''create table if not exists {outtbl} (
    r integer,
    m integer,
    ix integer,
    iy integer,

    m2 integer,
    x float,
    y float,

    dx float,
    dy float,
    sx float,
    sy float,
    snr float,

    dxb float,
    dyb float,
    sxb float,
    syb float,
    snrb float )''')

def transrunmont(r, m, ix, iy):
    print(f'Working on r{r} m{m} ix{ix},iy{iy}')
    img = rawimage.partialq5img(r, m, 0, ix, iy)
    x0,y0 = db.sel(f'''select x,y from {roughtbl}
    where r={r} and m={m} and s=0''')[0]
    s1 = ri.nslices(r-1) - 1
    if r-1==35:
        s1 -= 1
    mm1,xx1,yy1 = db.vsel(f'''select m,x,y from {roughtbl} 
    where r={r-1} and s={s1}
    order by m''')
    dx,dy = db.sel(f'select dx+dxb, dy+dyb from {intbl} where r2={r}')[0]
    dx *= inq/outq
    dy *= inq/outq

    # Center of tile in run space
    xc = x0 + X*(ix+.5)
    yc = y0 + Y*(iy+.5)

    # Center of other tiles
    xx1c = xx1 + (X*2.5)
    yy1c = yy1 + (X*2.5)

    # Presumed best montage to match
    m1 = mm1[np.argmin((xc+dx - xx1c)**2 + (yc-dy - yy1c)**2)]
    # Position in that montage that should match
    x1c = xc-dx - xx1[m1]
    y1c = yc-dy - yy1[m1]
    ix1 = int(x1c/X+.5)
    iy1 = int(y1c/Y+.5)
    if ix1<=0:
        ix1=1
    elif ix1>=5:
        ix1=4
    if iy1<=0:
        iy1=1
    elif iy1>=5:
        iy1=4
    s = ri.nslices(r-1)-1
    if r-1 == 35:
      s -= 1
    img1 = rawimage.q5subimg2x2(r-1, m1, s, ix1, iy1)

    # Position in second image that should match center of first
    x1croi = x1c - (ix1-1)*X
    y1croi = y1c - (iy1-1)*Y

    if x1croi<0 or x1croi>2*X or y1croi<0 or y1croi>2*Y:
        # Have nothing to connect to
        db.exe(f'''insert into {outtbl}
            (r, m, ix, iy,
            m2, x, y,
            dx, dy, sx, sy, snr,
            dxb, dyb, sxb, syb, snrb)
            values
            ({r}, {m}, {ix}, {iy},
            {m1}, {x1c}, {y1c},
            0,0,0,0,0,
            0,0,0,0,0)''')
        return

    if PICTURES:
        qp.figure('/tmp/s1', 4, 4)
        qp.imsc(img, xx=np.arange(0,X), yy=-np.arange(0,Y))
        qp.pen('r', 2)
        qp.marker('+')
        qp.mark(X/2, -Y/2)
        qp.figure('/tmp/s2', 8, 8)
        qp.imsc(img1, xx=np.arange(0,X*2), yy=-np.arange(0,Y*2))
        qp.pen('r', 2)
        qp.marker('+')
        qp.mark(x1croi, -y1croi)

    SIZ = (512, 512)
    win1 = swiftir.extractStraightWindow(img, (X/2,Y/2), SIZ)
    apo1 = swiftir.apodize(win1)
    win2 = swiftir.extractStraightWindow(img1, (x1croi,y1croi), SIZ)
    apo2 = swiftir.apodize(win2)

    if PICTURES:
        qp.figure('/tmp/s1b', 4, 4)
        qp.imsc(apo1)
        
        qp.figure('/tmp/s2b', 4, 4)
        qp.imsc(apo2)

    (dx, dy, sx, sy, snr) = swiftir.swim(apo1, apo2)
    
    win2 = swiftir.extractStraightWindow(img1, (x1croi+dx, y1croi+dy), SIZ)
    apo2 = swiftir.apodize(win2)

    if PICTURES:
        qp.figure('/tmp/s1a', 4, 4)
        qp.imsc(apo1)
        
        qp.figure('/tmp/s2a', 4, 4)
        qp.imsc(apo2)
    
    (dxb, dyb, sxb, syb, snrb) = swiftir.swim(apo1, apo2)

    if dxb**2 + dyb**2 > 1:
        win2 = swiftir.extractStraightWindow(img1,
                                             (x1croi+dx+dxb, y1croi+dy+dyb),
                                             SIZ)
        apo2 = swiftir.apodize(win2)
        (dxc, dyc, sxc, syc, snrc) = swiftir.swim(apo1, apo2)
        if snrc>snrb:
            dx += dxb
            dy += dyb
            sx = sxb
            sy = syb
            snr = snrb
            dxb, dyb, sxb, syb, snrb = dxc, dyc, sxc, syc, snrc

    db.exe(f'''insert into {outtbl}
    (r, m, ix, iy,
    m2, x, y,
    dx, dy, sx, sy, snr,
    dxb, dyb, sxb, syb, snrb)
    values
    ({r}, {m}, {ix}, {iy},
    {m1}, {x1c}, {y1c},
    {dx}, {dy}, {sx}, {sy}, {snr}, 
    {dxb},{dyb},{sxb},{syb},{snrb})''')

def transrunmany(r, m):
    for ix in range(5):
        for iy in range(5):
            cnt = db.sel(f'''select count(*) from {outtbl}
            where r={r} and m={m} and ix={ix} and iy={iy}''')[0][0]
            if cnt==0:
                transrunmont(r, m, ix, iy)

maketable()
                
fac = factory.Factory(nthreads)
for r0 in range(1, ri.nruns()):
    r = r0+1
    cnt = db.sel(f'select count(*) from {outtbl} where r={r}')[0][0]
    if cnt<ri.nmontages(r)*IX*IY:
        for m in range(ri.nmontages(r)):
            cnt = db.sel(f'''select count(*) from {outtbl}
            where r={r} and m={m}''')[0][0]
            if cnt < IX*IY:
                fac.request(transrunmany, r, m)
fac.shutdown()

    
