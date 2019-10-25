#!/usr/bin/python3

# I will create a db table with (r, m, s, ix, iy),
# (dx, dy, sx, sy, snr), and (dxb, dyb, sxb, syb, snrb).
# The logic will be that s pixel (x,y) fits to slice (r, m, s-1)
# pixel (x+dx+dxb,y+dy+dyb).
# It will be understood that the "b" alignment was centered at
# (X/2-dx/2, Y/2-dy/2) of the "s" image.
# I'll call the table RELMONTALIGNQ5.
# I will do the calculation for each of the 25 subtiles in Q5.
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

nthreads = 12

db = aligndb.DB()

def droptable():
    db.exe('''drop table relmontalignq5''')

def maketable():
    db.exe('''create table if not exists relmontalignq5 (
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
    if neighborhoodimg is None:
        db.exe(f'''insert into relmontalignq5 
        (r,m,s,ix,iy,
        dx,dy,sx,sy,snr, dxb,dyb,sxb,syb,snrb)
        values
        ({r},{m},{s},{ix},{iy},
        0,0,0,0,1000,
        0,0,0,0,1000)''')
        return
    
    Y,X = tileimg.shape
    SIZ = (512, 512)
    win1 = swiftir.extractStraightWindow(tileimg, (X/2,Y/2), SIZ)
    win2 = swiftir.extractStraightWindow(neighborhoodimg, (X/2,Y/2), SIZ)

    apo1 = swiftir.apodize(win1)
    apo2 = swiftir.apodize(win2)
    (dx, dy, sx, sy, snr) = swiftir.swim(apo1, apo2)

    win1 = swiftir.extractStraightWindow(tileimg, (X/2-dx/2,Y/2-dy/2), SIZ)
    win2 = swiftir.extractStraightWindow(neighborhoodimg, (X/2+dx/2,Y/2+dy/2), SIZ)

    apo1 = swiftir.apodize(win1)
    apo2 = swiftir.apodize(win2)
    (dxb, dyb, sxb, syb, snrb) = swiftir.swim(apo1, apo2)

    db.exe(f'''insert into relmontalignq5 
    (r,m,s,ix,iy,
    dx,dy,sx,sy,snr, dxb,dyb,sxb,syb,snrb)
    values
    ({r},{m},{s},{ix},{iy},
    {dx},{dy},{sx},{sy},{snr}, 
    {dxb},{dyb},{sxb},{syb},{snrb})''')
        
def alignmanysubtiles(r, m, ix, iy):
    cnt = db.sel(f'''select count(1) from relmontalignq5 
    where r={r} and m={m} and ix={ix} and iy={iy}''')[0][0]
    if cnt==ri.nslices(r):
        return
    def loader(subtileid):
        r,m,s,ix,iy = subtileid
        try:
            img = rawimage.partialq5img(r,m,s,ix,iy)
        except Exception as e:
            print(e)
            img = np.zeros((684,684), dtype=np.uint8) + 128
        return img 
            
    def saver(subtileid, tileimg, neighborhoodimg):
        r,m,s,ix,iy = subtileid
        alignsubtiles(r, m, s, ix, iy, tileimg, neighborhoodimg)

    ssdone = set()
    for row in db.sel(f'''select s from relmontattouchq5
           where r={r} and m={m} and ix={ix} and iy={iy}'''):
        ssdone.add(int(row[0]))

    lastimg = None
    for s in range(ri.nslices(r)):
        if s in ssdone:
            lastimg = None
        else:
            subtileid = (r,m,s,ix,iy)
            img = loader(subtileid)
            if s>0 and lastimg is None:
                lastimg = loader((r,m,s-1,ix,iy))
            saver(subtileid, img, lastimg)
            lastimg = img

fac = factory.Factory(nthreads)
    
def queuealignmanysubtiles(r, m, ix, iy):
    cnt = db.sel(f'''select count(1) from relmontalignq5 
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

print(f'Waiting for factory to complete tasks')
fac.shutdown()    
        
