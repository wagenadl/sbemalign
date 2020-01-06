#!/usr/bin/python3

import aligndb
import time
import sys
import traceback
import config
import rawimage
import pyqplot as qp
import numpy as np
import swiftir
import factory

outtbl = 'interrunq25'

nthreads = 6
PICTURES = False

db = aligndb.DB()
ri = db.runinfo()

db.exe(f'''create table if not exists {outtbl} (
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

def interrun(r2):
    r1 = r2 - 1
    n = db.sel(f'select count(*) from {outtbl} where r1={r1} and r2={r2}')
    if n[0][0]>0:
        return
    print(f'Working on R{r1}:{r2}')
    s1 = ri.nslices(r1)-1
    if r1==35:
      s1 -= 1 
    img1 = rawimage.loadimage(f'{config.root}/slicesq25/R{r1}S{s1}.jpg')
    if img1 is None:
        raise Exception(f'Failed to read last image for R{r1}')
        
    img2 = rawimage.loadimage(f'{config.root}/slicesq25/R{r2}S0.jpg')
    if img2 is None:
        raise Exception(f'Failed to read first image for R{r2}')

    Y1,X1 = img1.shape
    Y2,X2 = img2.shape
    SIZ = (512,1024)
    
    win1 = swiftir.extractStraightWindow(img1, (X1/2,Y1/2), SIZ)
    win2 = swiftir.extractStraightWindow(img2, (X2/2,Y2/2), SIZ)
    apo1 = swiftir.apodize(win1)
    apo2 = swiftir.apodize(win2)

    if PICTURES:
        qp.figure('/tmp/s1', 4, 8)
        qp.imsc(apo1)
        
        qp.figure('/tmp/s2', 4, 8)
        qp.imsc(apo2)
    
    (dx, dy, sx, sy, snr) = swiftir.swim(apo1, apo2)
    
    win1 = swiftir.extractStraightWindow(img1, (X1/2-dx/2,Y1/2-dy/2), SIZ)
    win2 = swiftir.extractStraightWindow(img2, (X2/2+dx/2,Y2/2+dy/2), SIZ)
    apo1 = swiftir.apodize(win1)
    apo2 = swiftir.apodize(win2)

    if PICTURES:
        qp.figure('/tmp/s1a', 4, 8)
        qp.imsc(apo1)
        
        qp.figure('/tmp/s2a', 4, 8)
        qp.imsc(apo2)
    
    (dxb, dyb, sxb, syb, snrb) = swiftir.swim(apo1, apo2)

    M1 = ri.nmontages(r1)
    M2 = ri.nmontages(r2)
    print(f'R{r1}:{r2} [{M1}/{M2}] {dx:6.2f} {dy:6.2f} {snr:6.1f} | {dxb:6.2f} {dyb:6.2f} {snrb:6.1f}')

    db.exe(f'''insert into {outtbl}
      (r1,r2, dx,dy,sx,sy,snr, dxb,dyb,sxb,syb,snrb)
      values
      ({r1},{r2},
      {dx},{dy},{sx},{sy},{snr}, 
      {dxb},{dyb},{sxb},{syb},{snrb})''')


fac = factory.Factory(nthreads)
for r0 in range(1, ri.nruns()):
    r2 = r0+1
    fac.request(interrun, r2)
fac.shutdown()

    
