#!/usr/bin/python3

# I will create a db table with (r, m, s, ix, iy),
# (dx, dy, sx, sy, snr), and (dxb, dyb, sxb, syb, snrb).
# The logic will be that s pixel (x,y) fits to slice (r, m, s-1)
# pixel (x+dx+dxb,y+dy+dyb).
# It will be understood that the "b" alignment was centered at
# (X/2-dx/2, Y/2-dy/2) of the "s" image.
# I'll call the table RELMONTALIGNQ1.
# I will do the calculation for each of 15x15 1024x1024-pixel
# window at Q1.
# Windows are _centered_ at ix*512 in the unaligned imagery,
# i.e., at ix*512 + 102 in the 16-bit raw tiffs.
# We use ix =  3,  5,  7,  9, 11, 13, 15,
#             18, 20, 22, 24, 26, 28, 30.
# (The center is undersampled, because distortions are
# thought to be extremely small there.)
# Note that the corners and edges will be sampled by
# the interslice process, to be written.

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

ixx = [3,  5,  7,  9, 11, 13, 15,
       18, 20, 22, 24, 26, 28, 30 ]

db = aligndb.DB()

def droptable():
    db.exe('''drop table relmontalignq1''')

def maketable():
    db.exe('''create table if not exists relmontalignq1 (
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
        db.exe(f'''insert into relmontalignq1 
        (r,m,s,ix,iy,
        dx,dy,sx,sy,snr, dxb,dyb,sxb,syb,snrb)
        values
        ({r},{m},{s},{ix},{iy},
        0,0,0,0,1000,
        0,0,0,0,1000)''')
        return
    
    Y,X = tileimg.shape
    SIZ = (1024, 1024)
    win1 = swiftir.extractStraightWindow(tileimg, (X/2,Y/2), SIZ)
    win2 = swiftir.extractStraightWindow(neighborhoodimg, (X/2,Y/2), SIZ)

    apo1 = swiftir.apodize(win1)
    apo2 = swiftir.apodize(win2)
    (dx, dy, sx, sy, snr) = swiftir.swim(apo1, apo2)

    SIZ = (768, 768)
    win1 = swiftir.extractStraightWindow(tileimg, (X/2-dx/2,Y/2-dy/2), SIZ)
    win2 = swiftir.extractStraightWindow(neighborhoodimg, (X/2+dx/2,Y/2+dy/2), SIZ)

    apo1 = swiftir.apodize(win1)
    apo2 = swiftir.apodize(win2)
    (dxb, dyb, sxb, syb, snrb) = swiftir.swim(apo1, apo2)

    db.exe(f'''insert into relmontalignq1 
    (r,m,s,ix,iy,
    dx,dy,sx,sy,snr, dxb,dyb,sxb,syb,snrb)
    values
    ({r},{m},{s},{ix},{iy},
    {dx},{dy},{sx},{sy},{snr}, 
    {dxb},{dyb},{sxb},{syb},{snrb})''')
        
def alignmanysubtiles(r, m, ix, iy):
    cnt = db.sel(f'''select count(1) from relmontalignq1 
    where r={r} and m={m} and ix={ix} and iy={iy}''')[0][0]
    if cnt==ri.nslices(r):
        return
    def loader(subtileid):
        r,m,s,ix,iy = subtileid
        try:
            img = rawimage.q1roi(r,m,s,ix*512-512, iy*512-512, 1024, 1024)
        except Exception as e:
            print(e)
            img = np.zeros((1024, 1024), dtype=np.uint8) + 128
        return img 
            
    def saver(subtileid, tileimg, neighborhoodimg):
        r,m,s,ix,iy = subtileid
        alignsubtiles(r, m, s, ix, iy, tileimg, neighborhoodimg)
    if cnt>0:
        db.exe(f'''delete from relmontalignq1 
        where r={r} and m={m} and ix={ix} and iy={iy}''')
    lastimg = None
    for s in range(ri.nslices(r)):
        subtileid = (r,m,s,ix,iy)
        img = loader(subtileid)
        saver(subtileid, img, lastimg)
        lastimg = img

fac = factory.Factory(nthreads)
    
def queuealignmanysubtiles(r, m, ix, iy):
    cnt = db.sel(f'''select count(1) from relmontalignq1 
    where r={r} and m={m} and ix={ix} and iy={iy}''')[0][0]
    if cnt==ri.nslices(r):
        return
    fac.request(alignmanysubtiles, r, m, ix, iy)

def queuealignmontage(r, m):
    for ix in ixx:
        for iy in ixx:
            queuealignmanysubtiles(r, m, ix, iy)

maketable()
    
for r0 in range(ri.nruns()):
    r  = r0 + 1
    for m in range(ri.nmontages(r)):
        queuealignmontage(r, m)

print(f'Waiting for factory to complete tasks')
fac.shutdown()    
        
