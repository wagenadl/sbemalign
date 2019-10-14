#!/usr/bin/python3

# I will create a db table with (r, m, s),
# (dx, dy, sx, sy, snr), and (dxb, dyb, sxb, syb, snrb).
# The logic will be that s pixel (x,y) fits to center-of-run
# pixel(x+dx+dxb,y+dy+dyb).
# It will be understood that the "b" alignment was centered at
# (X/2-dx/2, Y/2-dy/2) of the "m1" image.
# I'll call the table MONTAGEALIGNQ25CO.
# This version works from the center out. It uses swiftir's buildout
# to reduce accumulation of errors.
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

nthreads = 10

db = aligndb.DB()

def droptable():
    db.exe('''drop table montagealignq25co''')

def maketable():
    db.exe('''create table if not exists montagealignq25co (
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
    if neighborhoodimg is None:
        db.exe(f'''insert into montagealignq25co 
        (r,m,s,
        dx,dy,sx,sy,snr, dxb,dyb,sxb,syb,snrb)
        values
        ({r},{m},{s},
        0,0,0,0,100,
        0,0,0,0,100)''')
        return tileimg
    
    Y,X = tileimg.shape
    SIZ = (512,512)
    win1 = swiftir.extractStraightWindow(tileimg, (X/2,Y/2), SIZ)
    win2 = swiftir.extractStraightWindow(neighborhoodimg, (X/2,Y/2), SIZ)
    print('win1: %.1f - %.1f. win2: %.1f - %.1f' %
          (np.min(win1), np.max(win1),
          np.min(win2), np.max(win2)))
    print(win1.dtype)
    print(win2.dtype)
    apo1 = swiftir.apodize(win1)
    apo2 = swiftir.apodize(win2)
    (dx, dy, sx, sy, snr) = swiftir.swim(apo1, apo2)

    win1 = swiftir.extractStraightWindow(tileimg, (X/2-dx,Y/2-dy), SIZ)
    apo1 = swiftir.apodize(win1)

    if False:
        qp.figure('/tmp/s1', 6, 3)
        qp.subplot(1,2,1)
        qp.imsc(apo1)
        qp.shrink()
        qp.subplot(1,2,2)
        qp.imsc(apo2)
        time.sleep(.25)
    
    (dxb, dyb, sxb, syb, snrb) = swiftir.swim(apo1, apo2)

    db.exe(f'''insert into montagealignq25co 
    (r,m,s,
    dx,dy,sx,sy,snr, dxb,dyb,sxb,syb,snrb)
    values
    ({r},{m},{s},
    {dx},{dy},{sx},{sy},{snr}, 
    {dxb},{dyb},{sxb},{syb},{snrb})''')

    dx += dxb
    dy += dyb
    return swiftir.extractStraightWindow(tileimg, (X/2-dx, Y/2-dy), (Y,X))
        
def alignmontage(r, m):
    def loader(tileid):
        r,m,s = tileid
        try:
            img = rawimage.q25img(r,m,s)
        except Exception as e:
            print(e)
            img = np.zeros((684,684), dtype=np.uint8) + 128
        return img 
            
    def saver(tileid, tileimg, neighborhoodimg, aux):
        r,m,s = tileid
        print(f'Working on r{r} m{m} s{s}')
        return aligntiles(r, m, s, tileimg, neighborhoodimg)
    db.exe(f'''delete from montagealignq25co where r={r} and m={m}''')
    ifns = []
    for s in range(ri.nslices(r)):
        tileid = (r,m,s)
        ifns.append(tileid)
    swiftir.buildout(ifns, loader=loader, transformer=saver, nbase=11)

fac = factory.Factory(nthreads)
    
def queuealignmontage(r, m):
    cnt = db.sel(f'''select count(1) from montagealignq25co 
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
        
