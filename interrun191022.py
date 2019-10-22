#!/usr/bin/python3

import aligndb
import time
import sys
import traceback

import rawimage
import warpq5run
import warp
import pyqplot as qp
import numpy as np
import swiftir

db = aligndb.DB()
ri = db.runinfo()

db.exe('''create table if not exists interrunq5 (
    r1 integer,
    r2 integer,

    dx float,
    dy float,
    sx float,
    sy float,
    snr float,

    dxb float,
    dyb float,
    sxb float,
    syb float,
    snrb float )''')

for r0 in range(ri.nruns()-1):
    r1 = r0+1
    r2 = r0+2
    n = db.sel(f'select count(*) from interrunq5 where r1={r1} and r2={r2}')
    if n[0][0]>0:
        continue
    print(f'Working on R{r1}:{r2}')
    
    img1 = warpq5run.warpedq5img(r1, ri.nslices(r1)-1)
    if img1 is None:
        print(f'Failed to read image for R{r1}')
        continue
    img2 = warpq5run.warpedq5img(r2, 0)
    if img2 is None:
        print(f'Failed to read image for R{r2}')
        continue

    Y1,X1 = img1.shape
    Y2,X2 = img2.shape
    SIZ = (2048,4096)
    
    win1 = swiftir.extractStraightWindow(img1, (X1/2,Y1/2), SIZ)
    win2 = swiftir.extractStraightWindow(img2, (X2/2,Y2/2), SIZ)
    apo1 = swiftir.apodize(win1)
    apo2 = swiftir.apodize(win2)
    
    qp.figure('/tmp/s1', 4, 8)
    qp.imsc(apo1)
    
    qp.figure('/tmp/s2', 4, 8)
    qp.imsc(apo2)
    
    (dx, dy, sx, sy, snr) = swiftir.swim(apo1, apo2)
    
    win1 = swiftir.extractStraightWindow(img1, (X1/2-dx/2,Y1/2-dy/2), SIZ)
    win2 = swiftir.extractStraightWindow(img2, (X2/2+dx/2,Y2/2+dy/2), SIZ)
    apo1 = swiftir.apodize(win1)
    apo2 = swiftir.apodize(win2)
    
    qp.figure('/tmp/s1a', 4, 8)
    qp.imsc(apo1)
    
    qp.figure('/tmp/s2a', 4, 8)
    qp.imsc(apo2)
    
    (dxb, dyb, sxb, syb, snrb) = swiftir.swim(apo1, apo2)

    M1 = ri.nmontages(r1)
    M2 = ri.nmontages(r2)
    print(f'R{r1}:{r2} [{M1}/{M2}] {dx:6.2f} {dy:6.2f} {snr:6.1f} | {dxb:6.2f} {dyb:6.2f} {snrb:6.1f}')

    db.exe(f'''insert into interrunq5
      (r1,r2, dx,dy,sx,sy,snr, dxb,dyb,sxb,syb,snrb)
      values
      ({r1},{r2},
      {dx},{dy},{sx},{sy},{snr}, 
      {dxb},{dyb},{sxb},{syb},{snrb})''')

