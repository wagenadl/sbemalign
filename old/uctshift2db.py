#!/usr/bin/python3

'''This calculates the shift between the beta alignment and the uct dataset.
This DEPENDS on beta2db as well as reslice. Specifically:
  - We must have the RUNS and BETAPOS tables in the db;
  - We must have the Q25 scaled tiles;
  - We must have the (Q100) scaled uCT slices.
This FILLS the UCTSHIFT table in the db with (r,m,s) and (dx,dy,sx,sy,snr)
data. A given dx value means that the estimated true center position of a
tile is betapos.xc + uctshift.dx.
See testemuct.py for experimental code.'''

import psycopg2
import time
import sys
import traceback
import queue
import threading

import swiftir
import pyqplot as qp
import numpy as np

import rawimage

q = 25
quctxy = 100
quctz = 10
nthr = 8

def tilefn(r, m, s):
    return rawimage.scaledtile(r, m, s, q)

def uctfn(z):
    uctroot =  '/lsi2/dw/170428/resliced-uct'
    return '%s/xy%iz%i-%i.ppm' % (uctroot, quctxy, quctz, z)

def z2uct(z):
    return z//quctz

db = psycopg2.connect(database='align170428')

def exe(sql, args=None):
    '''EXE - Execute a single SQL statement in its own transaction'''
    with db:
        with db.cursor() as c:
            c.execute(sql, args)

def nofail(sql, args=None):
    '''NOFAIL - Like EXE, but catches exceptions'''
    with db:
        with db.cursor() as c:
            try:
                c.execute(sql, args)
            except Exception as e:
                print(e)
                
def sel(sql, args=None):
    '''SEL - Run a SQL statement in its own transaction and return all
             fetched rows.'''
    with db:
        with db.cursor() as c:
            c.execute(sql, args)
            return c.fetchall()

######################################################################

exe('''create table if not exists uctshift (
  r integer, m integer, s integer,
  dx float, dy float, sx float, sy float, snr float, iter integer )''')

ri = sel('select r,M,S,z0 from runs')
nruns = len(ri)
nmont = {}
nslice = {}
z0 = {}
for row in ri:
    r = row[0]
    nmont[r] = row[1]
    nslice[r] = row[2]
    z0[r] = row[3]

######################################################################
def tryshift(fft_em2, uct, xc,yc, Ruct):
    xcuct = xc / quctxy
    ycuct = yc / quctxy
    roi_uct = swiftir.extractROI(uct, (int(xcuct - Ruct/2),
                                       int(ycuct - Ruct/2),
                                       Ruct, Ruct))
    if roi_uct is None:
      print('roi_uct is none:', xc, yc)
      raise Exception(f'roi_uct is none: {xc},{yc} [%s]' % str(uct.shape))
   
    apo_uct = swiftir.apodize(roi_uct)
    fft_uct = swiftir.fft(apo_uct)
    
    (dx, dy, sx, sy, snr) = swiftir.swim(fft_uct, fft_em2)
    return (dx, dy, sx, sy, snr)
    
def aligntile(uct, r, m, s):
    em = swiftir.loadImage(tilefn(r, m, s))
    if em is None:
        raise Exception('Cannot load tile %s' % tilefn(r,m,s))
    xc, yc = sel('select xc, yc from betapos where r=%s and m=%s and s=%s',
                 (r, m, s))[0]
    Y,X = em.shape
    Rem = 512
    Ruct = int(Rem * q/quctxy)
    roi_em = swiftir.extractROI(em, ((X-Rem)//2, (Y-Rem)//2, Rem, Rem))
    apo_em = swiftir.apodize(roi_em)
    keep = np.hstack((np.arange(Ruct//2),
                      np.arange(Ruct//2)+Rem-Ruct//2))
    fft_em = swiftir.fft(apo_em)
    fft_em2 = fft_em[keep,:,:]
    fft_em2 = fft_em2[:,keep,:]
    cdx = 0
    cdy = 0
    dx = 0
    dy = 0
    first = True
    iter = 0
    print(r,m,s)
    while first or (dx*dx + dy*dy > 2 and iter<5):
        (dx,dy,sx,sy,snr) = tryshift(fft_em2, uct, xc+cdx,yc+cdy, Ruct)
        cdx -= quctxy*dx
        cdy -= quctxy*dy
        if np.abs(cdx)>5000 or np.abs(cdy)>5000:
            print(f'Warning: failed to converge at R{r} M{m} S{s}')
            return (0,0,0,0,0,0)
        csx = quctxy*sx
        csy = quctxy*sy
        csnr = snr
        first = False
        iter += 1

    return (cdx, cdy, csx, csy, snr, iter)
    
def doslice(r, s):
    z = z0[r] + s
    
    mi = sel('select m from uctshift where r=%s and s=%s', (r,s))
    mm = {}
    for row in mi:
        mm[row[0]] = 1

    zuct = z2uct(z)
    ufn = uctfn(zuct) 
    uct = swiftir.loadImage(ufn)
    if uct is None:
        raise Exception(f'uct is none at r={r} s={s} z={z}/{zuct} fn={ufn}')
    for m in range(nmont[r]):
        if m not in mm:
            (dx, dy, sx, sy, snr, iter) = aligntile(uct, r, m, s)
            exe('''insert into uctshift (r,m,s, dx,dy,sx,sy, snr,iter)
            values (%s,%s,%s, %s,%s,%s,%s,%s,%s)''',
                (r,m,s, float(dx),float(dy),float(sx),float(sy),
                 float(snr),iter))
            
qu = queue.Queue()
def worker():
    while True:
        (r,s) = qu.get()
        try:
            print('Working on run %i, slice %i' % (r, s))
            doslice(r, s)
        except Exception as e:
            ei = sys.exc_info()
            tb = ei[2]
            print(traceback.print_tb(tb))
            print(e)
            print('FAILED TO GET SHIFT AT', r, s)
        qu.task_done()

threads = []
for n in range(nthr):
    t = threading.Thread(target=worker)
    t.start()
    threads.append(t)
        
def dorun(r):
    for s in range(nslice[r]):
        ndone = sel('select count(*) from uctshift where r=%s and s=%s', (r,s))[0][0]
        if ndone<nmont[r]:
            print('queuing slice', r, s)
            qu.put((r, s))
    print('waiting for workers to finish')
    qu.join()
            
######################################################################
    
for r in z0.keys():
    ndone = sel('select count(*) from uctshift where r=%s', (r,))[0][0]
    if ndone<nmont[r]*nslice[r]:
        print('queuing run', r)
        dorun(r)
