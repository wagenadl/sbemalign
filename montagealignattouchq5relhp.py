#!/usr/bin/python3

# I will create a db table with (r, m, s, ix, iy, x,y),
# (dx, dy, sx, sy, snr), and (dxb, dyb, sxb, syb, snrb).
# The logic will be that s pixel (x',y') fits to neighborhood
# pixel(x'-dx-dxb,y'-dy-dyb).
# It will be understood that the "b" alignment was centered at
# (x, y) of the image.
# Although this information is slightly redundant, I think it will be
# convenient to have it all.
# I'll call the table MONTAGEALIGNATTOUCHQ5RELHP.
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

SHOW = False

nthreads = 12

if SHOW:
    nthreads = 1

db = aligndb.DB()

def droptable():
    db.exe('''drop table montagealignattouchq5relhp''')

def maketable():
 
    db.exe('''create table if not exists montagealignattouchq5relhp (
    r integer,
    m integer,
    s integer,
    ix integer,
    iy integer,
    x float, 
    y float,
    m2 integer,

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

def alignsubtiles(subtileid, row, tileimg, neighborhoodimg):
    r, m, s, ix, iy = subtileid
    x,y,dx0,dy0,m2 = row
    # dx0,dy0 are positions of original subtile in slicealignq5
    print(f'Working on r{r} m{m} s{s} {ix},{iy} ({dx0},{dy0})')
    Y,X = tileimg.shape
    if dx0>dy0:
        SIZ = (X//4, Y//2)
    elif dy0>dx0:
        SIZ = (X//2, Y//4)
    else:
        SIZ = (X//2, Y//2)
    print(f'SIZ is {SIZ}')

    if neighborhoodimg is None:
        db.exe(f'''insert into montagealignattouchq5relhp 
        (r,m,s,ix,iy,x,y,m2,
        dx,dy,sx,sy,snr, dxb,dyb,sxb,syb,snrb)
        values
        ({r},{m},{s},{ix},{iy},{x},{y},{m2},
        0,0,0,0,100,
        0,0,0,0,100)''')
        return

    dxbase,dybase = baseshifts[subtileid] # These are the baseshifts
    # from montagealignq5coa
    x += dxbase
    y += dybase
    # So now we are nominally back at the pixel position where
    # the original cross-montage match was made

    win1 = swiftir.extractStraightWindow(tileimg, (x,y), SIZ)
    print('shape', win1.shape)
    win2 = swiftir.extractStraightWindow(neighborhoodimg, (x,y), SIZ)
    apo1 = swiftir.apodize(win1)
    apo2 = swiftir.apodize(win2)

    if SHOW:
        qp.figure('/tmp/s1', 8, 4)
        qp.subplot(1,2,1)
        qp.imsc(apo1)
        Y1,X1 = win1.shape
        qp.at(X1/2,Y1/2)
        qp.pen('b')
        qp.text(f'{ix},{iy} {dx0}/{dy0}: {SIZ}')
        qp.shrink(1,1)
        qp.subplot(1,2,2)
        qp.imsc(apo2)
        qp.shrink(1,1)
        #time.sleep(.5)
    
    (dx, dy, sx, sy, snr) = swiftir.swim(apo1, apo2)

    win1 = swiftir.extractStraightWindow(tileimg, (x-dx/2,y-dy/2), SIZ)
    win2 = swiftir.extractStraightWindow(neighborhoodimg, (x+dx/2,y+dy/2), SIZ)
    apo1b = swiftir.apodize(win1)
    apo2b = swiftir.apodize(win2)
    (dxb, dyb, sxb, syb, snrb) = swiftir.swim(apo1b, apo2b)

    dx += dxbase
    dy += dybase

    db.exe(f'''insert into montagealignattouchq5relhp 
    (r,m,s,ix,iy,x,y,m2,
    dx,dy,sx,sy,snr, dxb,dyb,sxb,syb,snrb)
    values
    ({r},{m},{s},{ix},{iy},{x},{y},{m2},
    {dx},{dy},{sx},{sy},{snr}, 
    {dxb},{dyb},{sxb},{syb},{snrb})''')
        
def alignmanysubtiles(r, m, ix, iy):
    def loader(subtileid):
        r,m,s,ix,iy = subtileid
        if subtileid in baseshifts:
            dx,dy = baseshifts[subtileid]
        else:
            dx,dy = db.sel(f'''select dx+dxb,dy+dyb from montagealignq5relhp
            where r={r} and m={m} and s={s} and ix={ix} and iy={iy}''')[0]
            baseshifts[subtileid] = (dx,dy)
        print(f'loading r{r} m{m} s{s} ix{ix} iy{iy}: %.1f %.1f' % (dx,dy))
        try:
            img = rawimage.partialq5img(r,m,s,ix,iy)
            Y,X = img.shape
            img = swiftir.extractStraightWindow(img, (X/2-dx, Y/2-dy), (X,Y))
        except Exception as e:
            print(e)
            img = np.zeros((684,684), dtype=np.uint8) + 128
        return img 
            
    def saver(subtileid, tileimg, neighborhoodimg, aux):
        r,m,s,ix,iy = subtileid
        rows1 = db.sel(f'''select x1,y1,dx0,dy0,m2 from slicealignq5pos 
        where r={r} and m1={m} and s={s} and ix1={ix} and iy1={iy}''')
        for row in rows1:
            alignsubtiles(subtileid, row, tileimg, neighborhoodimg)
        rows2 = db.sel(f'''select x2,y2,dx0,dy0,m1 from slicealignq5pos 
        where r={r} and m2={m} and s={s} and ix2={ix} and iy2={iy}''')
        for row in rows2:
            alignsubtiles(subtileid, row, tileimg, neighborhoodimg)

    db.exe(f'''delete from montagealignattouchq5relhp 
    where r={r} and m={m} and ix={ix} and iy={iy}''')
    subtileids = []
    for s in range(ri.nslices(r)):
        subtileids.append((r,m,s,ix,iy))
    swiftir.buildout(subtileids, nbase=11,
                     loader=loader, saver=saver)

fac = factory.Factory(nthreads)
    
def queuealignmanysubtiles(r, m, ix, iy):
    cnt = db.sel(f'''select count(1) from montagealignattouchq5relhp 
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

#droptable()
maketable()
    
for r0 in range(ri.nruns()):
    r  = r0 + 1
    for m in range(ri.nmontages(r)):
        queuealignmontage(r, m)

print(f'Waiting for factory to end')
fac.shutdown()    
        
