#!/usr/bin/python3

# I will create a db table with (r, m, s, m2, ii, x,y),
# (dx, dy, sx, sy, snr), and (dxb, dyb, sxb, syb, snrb).
# The logic will be that s pixel (x',y') fits to neighborhood
# pixel(x'+dx+dxb,y'+dy+dyb). Neighborhood is s-1.
# It will be understood that the "b" alignment was centered at
# (x, y) of the image.
# Although this information is slightly redundant, I think it will be
# convenient to have it all.
# I'll call the table RELMONTATTOUCHQ5.
# I will do the calculation for each of the positions defined by the
# slicealignq5 table.
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
    db.exe('''drop table relmontattouchq5''')

def maketable():
 
    db.exe('''create table if not exists relmontattouchq5 (
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

def alignsubtiles(subtileid, m2, row, tileimg, neighborhoodimg):
    r, m, s, ii, ix, iy = subtileid
    x, y = row
    # dx0,dy0 are positions of original subtile in slicealignq5
    print(f'Working on r{r} m{m}:m2{m} s{s} {ii}')
    Y,X = tileimg.shape
    if m2==ri.mleft(r,m) or m2==ri.mright(r,m):
        SIZ = (X//4, Y//2)
    elif m2==ri.mabove(r,m) or m2==ri.mbelow(r,m):
        SIZ = (X//2, Y//4)
    else:
        raise ValueError('Surprising pair')
        SIZ = (X//2, Y//2)

    if neighborhoodimg is None:
        db.exe(f'''insert into relmontattouchq5 
        (r,m,s,ix,iy,x,y,m2,
        dx,dy,sx,sy,snr, dxb,dyb,sxb,syb,snrb)
        values
        ({r},{m},{s},{ix},{iy},{x},{y},{m2},
        0,0,0,0,0,
        0,0,0,0,0)''')
        return

    win1 = swiftir.extractStraightWindow(tileimg, (x,y), SIZ)
    # print('shape', win1.shape)
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

    db.exe(f'''insert into relmontattouchq5 
    (r,m,s,ix,iy,x,y,m2,
    dx,dy,sx,sy,snr, dxb,dyb,sxb,syb,snrb)
    values
    ({r},{m},{s},{ix},{iy},{x},{y},{m2},
    {dx},{dy},{sx},{sy},{snr}, 
    {dxb},{dyb},{sxb},{syb},{snrb})''')
        
def alignmanysubtiles(r, m, m2, ii, iifix):
    ix,iy = iifix
    if ix is None:
        ix = ii
    if iy is None:
        iy = ii
    def loader(subtileid):
        r,m,s,ii,ix,iy = subtileid
        print(f'loading r{r} m{m} s{s} ix{ix} iy{iy}')
        try:
            img = rawimage.partialq5img(r,m,s,ix,iy)
        except Exception as e:
            print(e)
            img = np.zeros((684,684), dtype=np.uint8) + 128
        return img 
            
    def saver(subtileid, tileimg, neighborhoodimg):
        r,m,s,ii,ix,iy = subtileid
        rows1 = db.sel(f'''select
          x1-dx/2-dxb/2-dxc/2,y1-dy/2-dyb/2-dyc/2
          from slicealignq5
          where r={r} and m1={m} and m2={m2} and s={s} and ii={ii}''')
        for row in rows1:
            alignsubtiles(subtileid, m2, row, tileimg, neighborhoodimg)
        rows2 = db.sel(f'''select
        x2+dx/2+dxb/2+dxc/2,y2+dy/2+dyb/2+dyc/2
        from slicealignq5pos 
        where r={r} and m2={m} and m2={m2} and s={s} and ii={ii}''')
        for row in rows2:
            alignsubtiles(subtileid, m2, row, tileimg, neighborhoodimg)

    # Figure out what slices have to be done
    ssdone = set()
    for row in db.sel(f'''select s from relmontattouchq5
       where r={r} and m={m} and m2={m2} and ii={ii}'''):
        s = row[0]
        m2 = row[1]
        ssdone.add(s)
    
    lastimg = None
    for s in range(ri.nslices(r)):
        if s in ssdone:
            lastimg = None
        else:
            img = loader((r,m,s,ii,ix,iy))
            if s>0 and lastimg is None:
                lastimg = loader((r,m,s-1,ii,ix,iy))
            saver((r,m,s,ii,ix,iy), img, lastimg)
            lastimg = img

fac = factory.Factory(nthreads)
    
def queuealignmanysubtiles(r, m, m2, ii, iifix):
    cnt = db.sel(f'''select count(1) from relmontattouchq5 
        where r={r} and m={m} and m2={m2} and ii={ii}''')[0][0]
    cnt1 = db.sel(f'''select count(1) from slicealignq5
        where r={r} and m1={m} and m2={m2} and ii={ii}''')[0][0]
    cnt2 = db.sel(f'''select count(1) from slicealignq5
        where r={r} and m2={m} and m1={m2} and ii={ii}''')[0][0]
    if cnt == cnt1 + cnt2:
        return
    fac.request(alignmanysubtiles, r, m, m2, ii, iifix)

def queuealignmontagepair(r, m, m2, iifix):
    cnt = db.sel(f'''select count(1) from {outtbl}
    where r={r} and m={m} and m2={m2}''')
    cnt1 = db.sel(f'''select count(1) from slicealignq5
    where r={r} and m1={m} and m2={m2}''')[0][0]
    cnt2 = db.sel(f'''select count(1) from slicealignq5
    where r={r} and m2={m} and m1={m2}''')[0][0]
    if cnt==cnt1+cnt2:
        return
    
    for ii in range(5):
        queuealignmanysubtiles(r, m, m2, ii, iifix)
    
    
def queuealignmontage(r, m):
    mleft = ri.mleft(r, m)
    mright = ri.mright(r, m)
    mabove = ri.mabove(r, m)
    mbelow = ri.mbelow(r, m)
    if mleft is not None:
        queuealignmontagepair(r, m, mleft, (0, None))
    if mright is not None:
        queuealignmontagepair(r, m, mright, (4, None))
    if mabove is not None:
        queuealignmontagepair(r, m, mabove, (None, 0))
    if mbelow is not None:
        queuealignmontagepair(r, m, mbelow, (None, 4))

#droptable()
maketable()
    
for r0 in range(ri.nruns()):
    r  = r0 + 1
    cnt = db.sel(f'''select count(1) from relmontattouchq5 
    where r={r}''')[0][0]
    cnt1 = db.sel(f'''select count(1) from slicealignq5
    where r={r}''')[0][0]
    if cnt==2*cnt1:
        continue
    for m in range(ri.nmontages(r)):
        queuealignmontage(r, m)

print(f'Waiting for factory to end')
fac.shutdown()    
        
