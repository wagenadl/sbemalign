
import aligndb
import time
import sys
import traceback

import swiftir
import pyqplot as qp
import numpy as np

import rawimage
import factory

nthreads = 1

db = aligndb.DB()
ri = db.runinfo()

ixx = [3,  5,  7,  9, 11, 13, 15,
       18, 20, 22, 24, 26, 28, 30 ]
R = 17100//5//5
MARG = (17100 - 33*512)//2

def droptables():
    db.nofail('''drop table crossalignq1''')
    db.nofail('''drop table edgealignq1''')

def maketables():
    db.exe('''create table if not exists crossalignq1 (
    r integer,
    m1 integer,
    m2 integer,
    s integer,
    ii integer,
    x1 float,
    y1 float,
    x2 float,
    y2 float,
    dx float,
    dy float,
    sx float,
    sy float,
    snr float,
    dxb float,
    dyb float,
    sxb float,
    syb float,
    snrb float
    )''')

    db.exe('''create table if not exists edgealignq1 (
    r integer,
    m integer,
    m2 integer,
    s integer,
    ii integer,
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
    snrb float
    )''')

def edgecoord(r, m1, m2):
    cc = 512*np.array(ixx)
    C = ri.ncolumns(r)
    if m2==m1+C:
        # Two montages top-to-bottom
        dx,y1,y2 = db.vsel(f'''select avg(x2-x1),avg(y1),avg(y2)
          from slicealignq5pos 
          where r={r} and m1={m1} and m2={m2}''')
        # Take average and upscale to Q1
        y1 = int(5 * y1 + 4*R) - MARG
        y2 = int(5 * y2) - MARG
        dx = int(5 * dx/2)
        return (cc-dx, y1+0*cc, cc+dx, y2+0*cc)
    elif m2==m1+1 and m2//C==m1//C:
        # Two montages left-to-right
        dy,x1,x2 = db.vsel(f'''select avg(y2-y1),avg(x1),avg(x2)
          from slicealignq5pos 
          where r={r} and m1={m1} and m2={m2}''')
        # Take average and upscale to Q1
        x1 = int(5 * x1 + 4*R) - MARG
        x2 = int(5 * x2) - MARG
        dy = int(5 * dy/2)
        return (x1+0*cc, cc-dy, x2+0*cc, cc+dy)
    else:
        raise ValueError(f'Coord request for nonexistent edge: R{r} M{m1}:{m2}')

def loadpatch(r, m, s, x, y):
    return rawimage.q1roi(r, m, s, x-512, y-512, 1024, 1024)

def alignone(img1, img2):
    Y,X = img1.shape        
    apo1 = swiftir.apodize(img1)
    apo2 = swiftir.apodize(img2)
    (dx, dy, sx, sy, snr) = swiftir.swim(apo1, apo2)
    SIZ = (768, 768)
    win1 = swiftir.extractStraightWindow(img1, (X/2-dx/2,Y/2-dy/2), SIZ)
    win2 = swiftir.extractStraightWindow(img2, (X/2+dx/2,Y/2+dy/2), SIZ)
    apo1b = swiftir.apodize(win1)
    apo2b = swiftir.apodize(win2)
    (dxb, dyb, sxb, syb, snrb) = swiftir.swim(apo1b, apo2b)
    return ((dx, dy, sx, sy, snr), (dxb, dyb, sxb, syb, snrb))

def aligncross(r, m1, m2, s, ii, x1, y1, x2, y2, img1, img2):
    ali, alib = alignone(img1, img2)
    dx, dy, sx, sy, snr = ali
    dxb, dyb, sxb, syb, snrb = alib
    sql = f'''insert into crossalignq1
             (r, m1, m2, s,
             ii, x1, y1, x2, y2,
             dx, dy, sx, sy, snr,
             dxb, dyb, sxb, syb, snrb)
             values
             ({r}, {m1}, {m2}, {s},
             {ii}, {x1}, {y1}, {x2}, {y2},
             {dx}, {dy}, {sx}, {sy}, {snr},
             {dxb}, {dyb}, {sxb}, {syb}, {snrb})'''
    return sql

def alignedge(r, m, m2, s, ii, x, y, img1, img2):
    ali, alib = alignone(img1, img2)
    dx, dy, sx, sy, snr = ali
    dxb, dyb, sxb, syb, snrb = alib
    sql = f'''insert into edgealignq1
             (r, m, m2, s,
             ii, x, y,
             dx, dy, sx, sy, snr,
             dxb, dyb, sxb, syb, snrb)
             values
             ({r}, {m1}, {m2}, {s},
             {ii}, {x}, {y},
             {dx}, {dy}, {sx}, {sy}, {snr},
             {dxb}, {dyb}, {sxb}, {syb}, {snrb})'''
    return sql

def alignsubtiles(r, m1, m2, ii, edgec):
    x1, y1, x2, y2 = edgec
    print(f'Aligning subtiles R{r} M{m1}:{m2} ({ii})')
    print('coords', x1, y1, x2, y2)
    img1 = None
    img2 = None
    for s in range(ri.nslices(r)):
        print(f'Working on R{r} M{m1}:{m2} S{s} ({ii})')
        if db.sel(f'''select count(*) from crossalignq1
           where r={r} and m1={m1} and m2={m2} and s={s} and ii={ii}''')[0][0]:
            img1 = None
            img2 = None
        else:
            sql = []
            img1a = img1
            img2a = img2
            img1 = loadpatch(r, m1, s, x1, y1)
            img2 = loadpatch(r, m2, s, x2, y2)
            sql.append(aligncross(r, m1, m2, s, ii, x1, y1, x2, y2, img1, img2))
            if s>0:
                if img1a is None:
                    img1a = loadpatch(r, m1, s-1, x1, y1)
                if img2a is None:
                    img2a = loadpatch(r, m2, s-1, x2, y2)
                sql.append(alignedge(r, m1, m2, s, ii, x1, y1, img1a, img1))
                sql.append(alignedge(r, m2, m1, s, ii, x2, y2, img2a, img2))
            with db.db:
                with db.db.cursor() as c:
                    for sq in sql:
                        c.execute(sq)
            
maketables()

fac = factory.Factory(nthreads)

def alignmontages(r, m1, m2):
    print(f'Considering run {r} M{m1}:{m2}')
    edgec = edgecoord(r, m1, m2)
    for k in range(len(ixx)):
        ii = ixx[k]
        cnt = db.sel(f'''select count(*) from crossalignq1
        where r={r} and m1={m1} and m2={m2} and ii={ii}''')[0][0]
        if cnt<ri.nslices(r):
            x1 = edgec[0][k]
            y1 = edgec[1][k]
            x2 = edgec[2][k]
            y2 = edgec[3][k]
            fac.request(alignsubtiles, r, m1, m2, ii, (x1,y1,x2,y2))

for r0 in range(ri.nruns()):
    r = r0 + 1
    print(f'Considering run {r}')
    C = ri.ncolumns(r)
    W = ri.nrows(r)
    for c in range(C):
        for w in range(W-1):
            # Gaps b/w top-bottom montages
            m1 = c + w*C
            m2 = c + (w+1)*C
            cnt = db.sel(f'''select count(*) from crossalignq1
               where r={r} and m1={m1} and m2={m2}''')[0][0]
            if cnt<len(ixx) * ri.nslices(r):
                alignmontages(r, m1, m2)

    for c in range(C-1):
        for w in range(W):
            # Gaps b/w left-right montages
            m1 = c + w*C
            m2 = c+1 + w*C
            cnt = db.sel(f'''select count(*) from crossalignq1
               where r={r} and m1={m1} and m2={m2}''')[0][0]
            if cnt<len(ixx) * ri.nslices(r):
                alignmontages(r, m1, m2)

print('waiting for factory to complete tasks')            
fac.shutdown()                        
