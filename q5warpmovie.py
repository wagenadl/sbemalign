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

try:
    os.mkdir('/tmp/movie')
except:
    pass

for r0 in range(ri.nruns()):
    r = r0 + 1
    try:
        os.mkdir(f'/tmp/movie/R{r}')
        os.mkdir(f'/tmp/movie/aR{r}')
    except Exception as e:
        print(e)
        continue

    print(f'Working on R{r}')
    
    S = ri.nslices(r)

    for s in range(S):
        print(f'Loading R{r} S{s}')
        img = warpq5run.warpedq5img(r, s)
        Y,X = img.shape
        win = swiftir.extractStraightWindow(img, (X/2,Y/2), (1024,1024))
        cv2.imwrite(f'/tmp/movie/R{r}/S{s}.jpg', win)
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
        cv2.imwrite(f'/tmp/movie/aR{r}/S{s}.jpg', win)
    os.system(f'ffmpeg -i "/tmp/movie/R{r}/S%d.jpg"  -c:v libx265 /lsi2/dw/170428/q5movies/R{r}.mp4')
    os.system(f'ffmpeg -i "/tmp/movie/aR{r}/S%d.jpg"  -c:v libx265 /lsi2/dw/170428/q5omovies/R{r}.mp4')
    os.system(f'rm /tmp/movie/R{r}/*.jpg')
    os.system(f'rm /tmp/movie/aR{r}/*.jpg')


