#!/usr/bin/python3

# I will create a db table with (r, m, s, ix, iy, x,y),
# (dx, dy, sx, sy, snr), and (dxb, dyb, sxb, syb, snrb).
# The logic will be that s pixel (x',y') fits to neighborhood
# pixel(x'-dx-dxb,y'-dy-dyb).
# It will be understood that the "b" alignment was centered at
# (x, y) of the image.
# Although this information is slightly redundant, I think it will be
# convenient to have it all.
# I'll call the table MONTAGEALIGNATTOUCHQ5.
# I will do the calculation for each of the positions defined by the
# slicealignq5pos table.
# Note that I am _not_ scaling the dx, dy, etc. coordinates back up to Q1.
# That's why "Q5" is in the table name.

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

def droptable():
    db.exe('''drop table montagealignattouchq5''')

def maketable():
    db.exe('''create table if not exists montagealignattouchq5 (
    r integer,
    m integer,
    s integer,
    ix integer,
    iy integer,
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

ri = db.runinfo()

def alignsubtiles(r, m, s, ix, iy, x, y, tileimg, neighborhoodimg):
    print(f'Working on r{r} m{m} s{s} {ix},{iy}')
    Y,X = tileimg.shape
    dx0,dy0,dxb0,dyb0 = db.sel(f'''select dx,dy,dxb,dyb from montagealignq5 
    where r={r} and m={m} and s={s} and ix={ix} and iy={iy}''')[0]
    dx0 += dxb0
    dy0 += dyb0
    SIZ = (X//4, Y//4)
    qp.figure('/tmp/s1', 8, 4)
    win1 = swiftir.extractStraightWindow(tileimg, (x,y), SIZ)
    win2 = swiftir.extractStraightWindow(neighborhoodimg, (x+dx0,y+dy0), SIZ)
    qp.subplot(1,2,1)
    qp.imsc(win1)
    qp.shrink()
    qp.subplot(1,2,2)
    qp.imsc(win2)
    qp.shrink()
    
    apo1 = swiftir.apodize(win1)
    apo2 = swiftir.apodize(win2)
    (dx, dy, sx, sy, snr) = swiftir.swim(apo1, apo2)
    #dx += dx0
    #dy += dy0

    win1 = swiftir.extractStraightWindow(tileimg, (x,y), SIZ)
    win2 = swiftir.extractStraightWindow(neighborhoodimg, (x+dx,y+dy), SIZ)
    apo1b = swiftir.apodize(win1)
    apo2b = swiftir.apodize(win2)
    (dxb, dyb, sxb, syb, snrb) = swiftir.swim(apo1b, apo2b)

    db.exe(f'''insert into montagealignattouchq5 
    (r,m,s,ix,iy,x,y,
    dx,dy,sx,sy,snr, dxb,dyb,sxb,syb,snrb)
    values
    ({r},{m},{s},{ix},{iy},{x},{y},
    {dx},{dy},{sx},{sy},{snr}, 
    {dxb},{dyb},{sxb},{syb},{snrb})''')
        
def alignmanysubtiles(r, m, ix, iy):
    def loader(subtileid):
        r,m,s,ix,iy = subtileid
        try:
            img = rawimage.partialq5img(r,m,s,ix,iy)
        except Exception as e:
            print(e)
            img = np.zeros((684,684), dtype=np.uint8) + 128
        return img 
            
    def saver(neighborhoodimg, subtileid, tileimg):
        r,m,s,ix,iy = subtileid
        rows1 = db.sel(f'''select x1,y1 from slicealignq5pos 
        where r={r} and m1={m} and s={s} and ix1={ix} and iy1={iy}''')
        for row in rows1:
            x,y = row
            alignsubtiles(r, m, s, ix, iy, x, y, tileimg, neighborhoodimg)
        rows2 = db.sel(f'''select x2,y2 from slicealignq5pos 
        where r={r} and m2={m} and s={s} and ix2={ix} and iy2={iy}''')
        for row in rows2:
            x,y = row
            alignsubtiles(r, m, s, ix, iy, x, y, tileimg, neighborhoodimg)

    db.exe(f'''delete from montagealignattouchq5 
    where r={r} and m={m} and ix={ix} and iy={iy}''')
    subtileids = []
    for s in range(ri.nslices(r)):
        subtileids.append((r,m,s,ix,iy))
    swiftir.remod(subtileids, halfwidth=10, topbot=True,
                  loader=loader, saver=saver)

fac = factory.Factory(nthreads)
    
def queuealignmanysubtiles(r, m, ix, iy):
    cnt = db.sel(f'''select count(1) from montagealignattouchq5 
    where r={r} and m={m} and ix={ix} and iy={iy}''')[0][0]
    cnt1 = db.sel(f'''select count(1) from slicealignq5
    where r={r} and m1={m} and ix1={ix} and iy1={iy}''')[0][0]
    cnt2 = db.sel(f'''select count(1) from slicealignq5
    where r={r} and m2={m} and ix2={ix} and iy2={iy}''')[0][0]
    if cnt==cnt1+cnt2:
        return
    fac.request(alignmanysubtiles, r, m, ix, iy)

def queuealignmontage(r, m):
    for ix in range(5):
        for iy in range(5):
            queuealignmanysubtiles(r, m, ix, iy)

droptable()
maketable()
    
for r0 in range(ri.nruns()):
    r  = r0 + 1
    for m in range(ri.nmontages(r)):
        queuealignmontage(r, m)

print(f'Waiting for factory to end')
fac.shutdown()    
        
