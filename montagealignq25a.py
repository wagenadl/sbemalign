#!/usr/bin/python3

# I will create a db table with (r, m, s),
# (dx, dy, sx, sy, snr), and (dxb, dyb, sxb, syb, snrb).
# The logic will be that s pixel (x,y) fits to center-of-run
# pixel(x+dx+dxb,y+dy+dyb).
# It will be understood that the "b" alignment was centered at
# (X/2-dx/2, Y/2-dy/2) of the "m1" image.
# Although this information is slightly redundant, I think it will be
# convenient to have it all.
# I'll call the table MONTAGEALIGNQ25A.
# Note that I am _not_ scaling the dx, dy, etc. coordinates back up to Q1.
# That's why "Q25" is in the table name.

import aligndb
import time
import sys
import traceback

import swiftir
import pyqplot as qp
import numpy as np

import rawimage
import factory

nthreads = 6

db = aligndb.DB()

def droptable():
    db.exe('''drop table montagealignq25a''')

def maketable():
    db.exe('''create table if not exists montagealignq25a (
    r integer,
    m integer,
    s integer,

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

def aligntiles(r, m, s, tileimg, neighborhoodimg):
    print(f'Working on r{r} m{m} s{s}')
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

    db.exe(f'''insert into montagealignq25a 
    (r,m,s,
    dx,dy,sx,sy,snr, dxb,dyb,sxb,syb,snrb)
    values
    ({r},{m},{s},
    {dx},{dy},{sx},{sy},{snr}, 
    {dxb},{dyb},{sxb},{syb},{snrb})''')
        
def alignmontage(r, m):
    def loader(tileid):
        r,m,s = tileid
        try:
            img = rawimage.q25img(r,m,s)
        except Exception as e:
            print(e)
            img = np.zeros((684,684), dtype=np.uint8) + 128
        return img 
            
    def saver(neighborhoodimg, tileid, tileimg):
        r,m,s = tileid
        aligntiles(r, m, s, tileimg, neighborhoodimg)
    db.exe(f'''delete from montagealignq25a where r={r} and m={m}''')
    S = ri.nslices(r)
    s0 = S//2
    img0 = rawimage.q25img(r,m,s0)
    for s in range(S):
        img = loader((r,m,s))
        saver(img0, (r,m,s), img)

fac = factory.Factory(nthreads)
    
def queuealignmontage(r, m):
    cnt = db.sel(f'''select count(1) from montagealignq25a 
    where r={r} and m={m}''')[0][0]
    if cnt==ri.nslices(r):
        return
    fac.request(alignmontage, r, m)

maketable()
    
for r0 in range(ri.nruns()):
    r  = r0 + 1
    for m in range(ri.nmontages(r)):
        queuealignmontage(r, m)

print(f'Waiting for factory to end')
fac.shutdown()    
        
