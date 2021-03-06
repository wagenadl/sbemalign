#!/usr/bin/python3

# I will create a db table with (r, m1, m2, s, ii),
# (x1,y1), (x2,y2), (dx, dy, sx, sy, snr), and (dxb, dyb, sxb, syb, snrb),
# and (dxc, dyc, sxc, syc, snrc).
# The logic will be that m2 pixel (x2+dx/2+dxb/2+dxc/2, y2+dy/2+dyb/2+dyc/2)
# fits to m1 pixel (x1-dx/2-dxb/2-dxc/2, y1-dy/2-dyb/2-dyc/2).
# The column ii counts where along the edge the alignment was done.
# Although this information is slightly redundant, I think it will be
# convenient to have it all. 
# I'll call the table SLICEALIGNQ5.
# I will do the calculation for each of the 5 edge subtiles in Q5.
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

nthreads = 6

db = aligndb.DB()

def droptable():
    db.exe('''drop table slicealignq5''')

def maketable():
    db.exe('''create table if not exists slicealignq5 (
    r integer,
    m1 integer,
    m2 integer,
    s integer,
    ii integer,
    x1 integer,
    y1 integer,
    x2 integer,
    y2 integer,
    dx float,
    dy float,
    sx float,
    sy float,
    snr float,
    dxb float,
    dyb float,
    sxb float,
    syb float,
    snrb float,
    dxc float,
    dyc float,
    sxc float,
    syc float,
    snrc float 
    )''')

ri = db.runinfo()

def alignsubtiles(r, m1, m2, s, ii, sidebyside):
    cnt = db.sel(f'''select count(1) from slicealignq5 
    where r={r} and m1={m1} and m2={m2} and s={s}
    and ii={ii}''')
    if cnt[0][0]>0:
        return
    print(f'Working on r{r} m{m1}:{m2} s{s} {ii}')
    if sidebyside:
        ix1 = 4
        iy1 = ii
        ix2 = 0
        iy2 = ii
    else:
        ix1 = ii
        iy1 = 4
        ix2 = ii
        iy2 = 0
    try:
        img1 = rawimage.partialq5img(r,m1,s,ix1,iy1)
    except Exception as e:
        print(e)
        img1 = np.zeros((684,684), dtype=np.uint8) + 128
    try:
        img2 = rawimage.partialq5img(r,m2,s,ix2,iy2)
    except Exception as e:
        print(e)
        img2 = np.zeros((684,684), dtype=np.uint8) + 128
    Y,X = img1.shape
    if sidebyside:
        if r==25 and s>140:
            dx0 = X//3
        else:
            dx0 = X//2
        x1l = X//2 + dx0//2 - X//4
        x2l = X//2 - dx0//2 - X//4
        x1r = x1l+X//2
        x2r = x2l+X//2
        img1 = img1[:,x1l:x1r]
        img2 = img2[:,x2l:x2r]
        dy0 = 0
        x1 = (x1l+x1r)//2 + X*4
        x2 = (x2l+x2r)//2
        y1 = y2 = Y//2 + Y*iy1
    else:
        img1 = img1[Y//2:,:]
        img2 = img2[:Y//2,:]
        dx0 = 0
        dy0 = Y//2
        y1t = Y//2 + dy0//2 - Y//4
        y2t = Y//2 - dy0//2 - Y//4
        y1b = y1t+Y//2
        y2b = y2t+Y//2
        x1 = x2 = X//2 + X*ix1
        y1 = (y1t+y1b)//2 + Y*4
        y2 = (y2t+y2b)//2
        #y1 = Y//2 - dy0/2 + Y*4
        #y2 = Y//2 + dy0/2
    Y,X = img1.shape        
    apo1 = swiftir.apodize(img1)
    apo2 = swiftir.apodize(img2)
    (dx, dy, sx, sy, snr) = swiftir.swim(apo1, apo2)

    win1 = swiftir.extractStraightWindow(img1, (X/2-dx/2,Y/2-dy/2),(X//2,Y//2))
    win2 = swiftir.extractStraightWindow(img2, (X/2+dx/2,Y/2+dy/2),(X//2,Y//2))
    apo1b = swiftir.apodize(win1)
    apo2b = swiftir.apodize(win2)
    (dxb, dyb, sxb, syb, snrb) = swiftir.swim(apo1b, apo2b)

    dx1 = dx+dxb
    dy1 = dy+dyb
    win1 = swiftir.extractStraightWindow(img1,(X/2-dx1/2,Y/2-dy1/2),(X//2,Y//2))
    win2 = swiftir.extractStraightWindow(img2,(X/2+dx1/2,Y/2+dy1/2),(X//2,Y//2))
    apo1b = swiftir.apodize(win1)
    apo2b = swiftir.apodize(win2)
    (dxc, dyc, sxc, syc, snrc) = swiftir.swim(apo1b, apo2b)

    #print(f'{r} {m1}:{m2} {s} {ix1},{iy1}:{ix2}:{iy2} {dx0}+{dx}+{dxb}+{dxc} {dy0}+{dy}+{dyb}+{dyc}')
    
    db.exe(f'''insert into slicealignq5 
    (r,m1,m2,s,ii, x1,y1, x2,y2,
    dx,dy,sx,sy,snr, dxb,dyb,sxb,syb,snrb, dxc,dyc,sxc,syc,snrc)
    values
    ({r},{m1},{m2},{s},{ii}, {x1},{y1}, {x2},{y2},
    {dx},{dy},{sx},{sy},{snr}, 
    {dxb},{dyb},{sxb},{syb},{snrb},
    {dxc},{dyc},{sxc},{syc},{snrc})''')

def alignsidebyside(r, m1, m2, s):
    for iy in range(5):
        try:
            alignsubtiles(r, m1, m2, s, iy, True)
        except Exception as e:
            print(e)
            print(f'Failed to align r{r} m{m1}:{m2} s{s} iy{iy}')

def alignabove(r, m1, m2, s):
    for ix in range(5):
        try:
            alignsubtiles(r, m1, m2, s, ix, False)
        except Exception as e:
            print(e)
            print(f'Failed to align r{r} m{m1}:{m2} s{s} ix{ix}')
        
def aligntiles(r, m1, m2, s):
    cnt = db.sel(f'''select count(1) from slicealignq5 
    where r={r} and m1={m1} and m2={m2} and s={s}''')
    if cnt[0][0]==5:
        return
    C = ri.ncolumns(r)
    if C==2 and (m1%2==0 and m2==m1+1):
        alignsidebyside(r, m1, m2, s)
    elif (C==1 and (m2==m1+1)) or (C==2 and (m2==m1+2)):
        alignabove(r, m1, m2, s)
    else:
        raise Exception(f'Cannot align MM {m1} and {m2} in R {r} S {s}')

def alignmanytiles(r, m1, m2):
    S = ri.nslices(r)
    for s in range(S):
        aligntiles(r, m1, m2, s)

fac = factory.Factory(nthreads)
        
def queuealignmanytiles(r, m1, m2):
    cnt = db.sel(f'''select count(1) from slicealignq5 
    where r={r} and m1={m1} and m2={m2}''')
    S = ri.nslices(r)
    if cnt[0][0]==5*S:
        return
    print(f'Requesting r{r} m{m1}:{m2}')
    fac.request(alignmanytiles, r, m1, m2)

maketable()
    
for r0 in range(ri.nruns()):
    r  = r0 + 1
    print(f'Considering run {r}')
    ROWS = ri.nrows(r)
    COLS = ri.ncolumns(r)
    cnt = db.sel(f'''select count(1) from slicealignq5 where r={r}''')[0][0]
    if COLS==1:
        if cnt<(ROWS-1)*5*ri.nslices(r):
            # Above
            for m in range(ROWS-1):
                queuealignmanytiles(r, m, m+1)
    elif COLS==2:
        if cnt<((ROWS-1)*COLS + ROWS*(COLS-1))*5*ri.nslices(r):
            # Above
            for row in range(ROWS-1):
                for col in range(COLS):
                    m = col + 2*row
                    queuealignmanytiles(r, m, m+2)
            # Side-by-side
            for row in range(ROWS):
                m = 2*row
                queuealignmanytiles(r, m, m+1)
    else:
        raise Exception(f'Not handling {COLS} columns')

print(f'Waiting for factory to end')
fac.shutdown()    
