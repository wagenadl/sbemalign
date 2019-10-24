#!/usr/bin/python3

import aligndb
import time
import sys
import traceback
import os
import rawimage
import warpq5run
import warp
import pyqplot as qp
import numpy as np
import swiftir
import cv2

db = aligndb.DB()
ri = db.runinfo()

def overviewfn(r):
    return f'/lsi2/dw/170428/q5omovies/R{r}.mp4'

def detailfn(r):
    return f'/lsi2/dw/170428/q5movies/R{r}.mp4'

def overviewtmp(r):
    return f'/tmp/movie/aR{r}'

def detailtmp(r):
    return f'/tmp/movie/R{r}'

def makedirs(r):
    if not os.path.exists('/tmp/movie'):
        os.mkdir('/tmp/movie')
    if not os.path.exists(detailtmp(r)):
        os.mkdir(detailtmp(r))
    if not os.path.exists(overviewtmp(r)):
        os.mkdir(overviewtmp(r))

for r0 in range(ri.nruns()):
    r = r0 + 1
    detfn = detailfn(r)
    ovfn = overviewfn(r)
    if os.path.exists(ovfn) and os.path.exists(detfn):
        continue

    print(f'Working on R{r}')
    
    S = ri.nslices(r)

    im0 = None
    dettmp = detailtmp(r)
    ovtmp = overviewtmp(r)

    makedirs(r)
    
    for s in range(S):
        print(f'Loading R{r} S{s}')
        img = warpq5run.warpedq5img(r, s)
        if img is None:
            print(f'Missing image for R{r} S{s} - blacking out')
            img = 0*im0
        Y,X = img.shape
        win = swiftir.extractStraightWindow(img, (X/2,Y/2), (1024,1024))
        cv2.imwrite(f'{dettmp}/S{s}.jpg', win)
        print(img.shape)
        SCL = 20
        NX = X//SCL
        NY = Y//SCL
        win = img[:SCL*NY, :SCL*NX]
        print(win.shape)
        win = np.reshape(win, (NY, SCL, NX, SCL))
        print(win.shape)
        win = np.mean(np.mean(win,3), 1)
        print(win.shape)
        win = np.transpose(win)
        cv2.imwrite(f'{ovtmp}/S{s}.jpg', win)
    os.system(f'ffmpeg -i "{dettmp}/S%d.jpg"  -c:v libx265 {detfn}')
    os.system(f'ffmpeg -i "{ovtmp}/S%d.jpg"  -c:v libx265 {ovfn}')
    os.system(f'rm -r {ovtmp}')
    os.system(f'rm -r {dettmp}')


