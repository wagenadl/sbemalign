#!/usr/bin/python3

# I will create a db table with (r, m, s, ix, iy),
# (dx, dy, sx, sy, snr), and (dxb, dyb, sxb, syb, snrb).
# The logic will be that s pixel (x,y) fits to neighborhood
# pixel(x+dx+dxb,y+dy+dyb).
# It will be understood that the "b" alignment was centered at
# (X/2-dx/2, Y/2-dy/2) of the "m" image.
# Although this information is slightly redundant, I think it will be
# convenient to have it all.
# I'll call the table MONTAGEALIGNQ5A.
# I will do the calculation for each of the 25 subtiles in Q5.
# Note that I am _not_ scaling the dx, dy, etc. coordinates back up to Q1.
# That's why "Q5" is in the table name.
# This is an improvement over MONTAGEALIGNQ5, in that it uses coarse alignment
# from MONTAGEALIGNQ25.

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
    db.exe('''drop table montagealignq5a''')

def maketable():
    db.exe('''create table if not exists montagealignq5a (
    r integer,
    m integer,
    s integer,
    ix integer,
    iy integer,

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

def alignsubtiles(r, m, s, ix, iy, tileimg, neighborhoodimg):
    print(f'Working on r{r} m{m} s{s} {ix},{iy}')
    Y,X = tileimg.shape

    apo1 = swiftir.apodize(tileimg)
    apo2 = swiftir.apodize(neighborhoodimg)
    (dx, dy, sx, sy, snr) = swiftir.swim(apo1, apo2)

    win1 = swiftir.extractStraightWindow(tileimg, (X/2-dx/2,Y/2-dy/2),
                                         (X*3//4,Y*3//4))
    win2 = swiftir.extractStraightWindow(neighborhoodimg, (X/2+dx/2,Y/2+dy/2),
                                         (X*3//4,Y*3//4))
    apo1b = swiftir.apodize(win1)
    apo2b = swiftir.apodize(win2)
    (dxb, dyb, sxb, syb, snrb) = swiftir.swim(apo1b, apo2b)

    db.exe(f'''insert into montagealignq5a 
    (r,m,s,ix,iy,
    dx,dy,sx,sy,snr, dxb,dyb,sxb,syb,snrb)
    values
    ({r},{m},{s},{ix},{iy},
    {dx},{dy},{sx},{sy},{snr}, 
    {dxb},{dyb},{sxb},{syb},{snrb})''')
        
def alignmanysubtiles(r, m, ix, iy):
    cnt = db.sel(f'''select count(1) from montagealignq5a
    where r={r} and m={m} and ix={ix} and iy={iy}''')[0][0]
    if cnt==ri.nslices(r):
        return
    def loader(subtileid):
        r,m,s,ix,iy = subtileid
        dx,dy,dxb,dyb = db.sel(f'''select dx,dy,dxb,dyb from montagealignq25a
        where r={r} and m={m} and s={s}''')[0]
        dx += dxb
        dy += dyb
        #dx *= 20/21
        #dy *= 20/21
        dx *= 5
        dy *= 5
        print(f'loading r{r} m{m} s{s} ix{ix} iy{iy}: %.1f %.1f' % (dx,dy))
        try:
            img = rawimage.partialq5img(r,m,s,ix,iy)
            Y,X = img.shape
            img = swiftir.extractStraightWindow(img, (X/2-dx, Y/2-dy), (X,Y))
            qp.figure('/tmp/s1', 3, 3)
            qp.imsc(img)
            qp.at(250,250)
            qp.pen('r')
            qp.text(f'r{r} m{m} s{s} ix{ix} iy{iy}')
            qp.text('dx = %.1f dy = %.1f' % (dx,dy), dy=12)
            time.sleep(.5)
        except Exception as e:
            print(e)
            img = np.zeros((684,684), dtype=np.uint8) + 128
        return img 
            
    def saver(neighborhoodimg, subtileid, tileimg):
        r,m,s,ix,iy = subtileid
        qp.figure('/tmp/s2', 6, 3)
        qp.subplot(1,2,1)
        qp.imsc(tileimg)
        qp.at(250,250)
        qp.pen('r')
        qp.text(f'r{r} m{m} s{s} ix{ix} iy{iy}')
        qp.shrink()
        qp.subplot(1,2,2)
        qp.imsc(neighborhoodimg)
        qp.shrink()
        time.sleep(.5)
        alignsubtiles(r, m, s, ix, iy, tileimg, neighborhoodimg)
    if cnt>0:
        db.exe(f'''delete from montagealignq5a 
        where r={r} and m={m} and ix={ix} and iy={iy}''')
    subtileids = []
    for s in range(ri.nslices(r)):
        subtileids.append((r,m,s,ix,iy))
    swiftir.remod(subtileids, halfwidth=10, topbot=True,
                  loader=loader, saver=saver)

fac = factory.Factory(nthreads)
    
def queuealignmanysubtiles(r, m, ix, iy):
    cnt = db.sel(f'''select count(1) from montagealignq5a
    where r={r} and m={m} and ix={ix} and iy={iy}''')[0][0]
    if cnt==ri.nslices(r):
        return
    fac.request(alignmanysubtiles, r, m, ix, iy)

def queuealignmontage(r, m):
    for ix in [2]: #range(5):
        for iy in [2]: #range(5):
            queuealignmanysubtiles(r, m, ix, iy)

droptable()
maketable()
    
for r0 in range(ri.nruns()):
    r  = r0 + 1
    for m in range(ri.nmontages(r)):
        queuealignmontage(r, m)

print(f'Waiting for factory to end')
fac.shutdown()    
        
