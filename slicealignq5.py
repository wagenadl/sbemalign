#!/usr/bin/python3

# I will create a db table with (r, m1, m2, s, ix1, iy1, ix2, iy2),
# (dx0, dy0), (dx, dy, sx, sy, snr), and (dxb, dyb, sxb, syb, snrb),
# and (dxc, dyc, sxc, syc, snrc).
# The logic will be that m2 pixel (x,y) fits to m1
# pixel(x+dx0-dx-dxb-dxc,y+dy0-dy-dyb-dyc).
# If dx0>0, it will be understood that the "b" alignment was centered at
# (dx0+X/4-dx/2, Y/2-dy/2) of the "m1" image, else if dy0>0, the "b" alignment
# is centered at (X/2-dx/2, dy0+Y/4-dy/2). It is always true that precisely
# one of dx0 or dy0 is nonzero.
# Although this information is slightly redundant, I think it will be
# convenient to have it all. 
# I'll call the table SLICEALIGNQ5.
# I will do the calculation for each of the 5 edge subtiles in Q5.
# Note that I am _not_ scaling the dx, dy, etc. coordinates back up to Q1.
# That's why "Q5" is in the table name.
# I should have picked less stupid sign conventions.
# Since in my implementation, dx0 is either 0 or X/2, and analogous for dy0,
# the annoying logic about can be written as
# (X/2+dx0/2-dx/2, Y/2+dy0/2-dy/2). This is an awful hack, but it does
# always work. Note that X = Y = 684.
# Similarly, the matching pixel in m2 is at (X/2-dx0/2+dx/2, Y/2-dy0/2+dy/2).
# Also worth noting, the matching region is X/2 × Y/4 pixels if dy0>0 or
# X/4 × Y/2 pixels if dx0>0. This can be written as
# X/2 - dx0/2 × Y/2 - dy0/2

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
    ix1 integer,
    ix2 integer,
    iy1 integer,
    iy2 integer,
    dx0 float,
    dy0 float,
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

def runinfo():
    ri = db.sel('select r,M,S,z0 from runs')
    class RI:
        def nruns(self):
            return self.R
        def nslices(self, r):
            return self.SS[r]
        def nmontages(self, r):
            return self.MM[r]
        def z0(self, r):
            return self.zz0[r]
        def ncolumns(self, r):
            if self.MM[r]<=3:
                return 1
            else:
                return 2
        def nrows(self, r):
            return self.nmontages(r) // self.ncolumns(r)
    res = RI()
    res.R = len(ri)
    res.MM = {}
    res.SS = {}
    res.zz0 = {}
    for row in ri:
        r = row[0]
        res.MM[r] = row[1]
        res.SS[r] = row[2]
        res.zz0[r] = row[3]
    return res
ri = runinfo()

def alignsubtiles(r, m1, m2, s, ix1, iy1, ix2, iy2, sidebyside):
    cnt = db.sel(f'''select count(1) from slicealignq5 
    where r={r} and m1={m1} and m2={m2} and s={s}
    and ix1={ix1} and iy1={iy1}''')
    if cnt[0][0]>0:
        return
    print(f'Working on {r} {m1}:{m2} {s} {ix1},{iy1}:{ix2},{iy2}')
    img1 = rawimage.partialq5img(r,m1,s,ix1,iy1)
    img2 = rawimage.partialq5img(r,m2,s,ix2,iy2)
    Y,X = img1.shape
    if sidebyside:
        img1 = img1[:,X//2:]
        img2 = img2[:,:X//2]
        dx0 = X//2
        dy0 = 0
    else:
        img1 = img1[Y//2:,:]
        img2 = img2[:Y//2,:]
        dx0 = 0
        dy0 = Y//2
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
    (r,m1,m2,s,ix1,iy1,ix2,iy2, dx0,dy0,
    dx,dy,sx,sy,snr, dxb,dyb,sxb,syb,snrb, dxc,dyc,sxc,syc,snrc)
    values
    ({r},{m1},{m2},{s},{ix1},{iy1},{ix2},{iy2}, {dx0},{dy0},
    {dx},{dy},{sx},{sy},{snr}, 
    {dxb},{dyb},{sxb},{syb},{snrb},
    {dxc},{dyc},{sxc},{syc},{snrc})''')

def alignsidebyside(r, m1, m2, s):
    for iy in range(5):
        alignsubtiles(r, m1, m2, s, 4,iy, 0,iy, True)

def alignabove(r, m1, m2, s):
    for ix in range(5):
        alignsubtiles(r, m1, m2, s, ix,4, ix,0, False)
        
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
    fac.request(alignmanytiles, r, m1, m2)

maketable()
    
for r0 in range(ri.nruns()):
    r  = r0 + 1
    ROWS = ri.nrows(r)
    COLS = ri.ncolumns(r)
    if COLS==1:
        # Above
        for m in range(ROWS-1):
            queuealignmanytiles(r, m, m+1)
    elif COLS==2:
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

