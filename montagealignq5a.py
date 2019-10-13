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
# Note that dx,dy are inclusive of the baseshift from MONTAGEALIGNQ5.

import aligndb
import time
import sys
import traceback

import swiftir
import pyqplot as qp
import numpy as np

import rawimage
import factory

nthreads = 12

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

baseshifts = {}

def alignsubtiles(r, m, s, ix, iy, tileimg, neighborhoodimg):
    print(f'Working on r{r} m{m} s{s} {ix},{iy}')
    Y,X = tileimg.shape
    SIZ = (X*3//4, Y*3//4)
    win1 = swiftir.extractStraightWindow(tileimg, (X/2,Y/2), SIZ)
    win2 = swiftir.extractStraightWindow(neighborhoodimg, (X/2,Y/2), SIZ)
    apo1 = swiftir.apodize(win1)
    apo2 = swiftir.apodize(win2)
    (dx, dy, sx, sy, snr) = swiftir.swim(apo1, apo2)

    win1 = swiftir.extractStraightWindow(tileimg, (X/2-dx/2,Y/2-dy/2), SIZ)
    win2 = swiftir.extractStraightWindow(neighborhoodimg, (X/2+dx/2,Y/2+dy/2),
                                         SIZ)
    apo1b = swiftir.apodize(win1)
    apo2b = swiftir.apodize(win2)
    (dxb, dyb, sxb, syb, snrb) = swiftir.swim(apo1b, apo2b)

    tileid = (r,m,s)
    dx0,dy0 = baseshifts[tileid]
    dx += dx0
    dy += dy0

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
        tileid = (r,m,s)
        if tileid in baseshifts:
            dx,dy = baseshifts[tileid]
        else:
            dx,dy = db.sel(f'''select dx+dxb,dy+dyb from montagealignq25
            where r={r} and m={m} and s={s}''')[0]
            dx *= 20/21
            dy *= 20/21
            dx *= 5
            dy *= 5
            baseshifts[tileid] = (dx,dy)
        print(f'loading r{r} m{m} s{s} ix{ix} iy{iy}: %.1f %.1f' % (dx,dy))
        try:
            img = rawimage.partialq5img(r,m,s,ix,iy)
            Y,X = img.shape
            img = swiftir.extractStraightWindow(img, (X/2-dx, Y/2-dy), (X,Y))
        except Exception as e:
            print(e)
            img = np.zeros((684,684), dtype=np.uint8) + 128
        return img 
            
    def saver(neighborhoodimg, subtileid, tileimg):
        r,m,s,ix,iy = subtileid
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
    for ix in range(5):
        for iy in range(5):
            queuealignmanysubtiles(r, m, ix, iy)

maketable()
    
for r0 in range(ri.nruns()):
    r  = r0 + 1
    for m in range(ri.nmontages(r)):
        queuealignmontage(r, m)

print(f'Waiting for factory to end')
fac.shutdown()    
        
