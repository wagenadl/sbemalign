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

r = 1
S = ri.nslices(r)
os.mkdir('/tmp/movie')

for s in range(S):
    print(f'Loading R{r} S{s}')
    img = warpq5run.warpedq5img(r, s)
    Y,X = img.shape
    win = swiftir.extractStraightWindow(img, (X/2,Y/2), (1024,1024))
    cv2.imwrite(f'/tmp/movie/R{r}S{s}.jpg', win)
    

