#!/usr/bin/python3

# This does a completely naive stacking of a montage with local
# averaged results

import rawimage
import aligndb
import swiftir
import pyqplot as qp
import time
import numpy as np
import skimage.filters

db = aligndb.DB()

r = 21
m = 4
q = 25


def loadtile(rms):
    # Loads a tile without any shift
    r,m,s = rms
    print('loadtile', r, m, s)
    return rawimage.loadimage(rawimage.scaledtile(r,m,s,q))

def swimtile(remodimg, rms, localimg):
    # REMODIMG is the leave-one-out average
    # LOCALIMG is the left-out image
    r,m,s = rms
    print('swimtile', r, m, s)
    qp.figure('/tmp/remod', 8, 8)
    qp.subplot(2,2,1)
    qp.imsc(remodimg)
    qp.subplot(2,3,4)
    qp.imsc(remodimg)
    
    Y,X = remodimg.shape
    R = 512
    marg = (X-R)//2
    remodimg = swiftir.extractStraightWindow(remodimg, siz=R)
    remodimg = swiftir.apodize(remodimg)
    sd = np.mean(skimage.filters.laplace(remodimg.astype('float'))**2)
    print(sd)
    db.exe('insert into studystackup0 (s,sd) values (%s,%s)', (s,float(sd)))
    localimg = swiftir.extractStraightWindow(localimg, siz=R)
    localimg = swiftir.apodize(localimg)
    qp.subplot(2,3,2)
    qp.imsc(remodimg)
    qp.subplot(2,3,5)
    qp.imsc(localimg)
    (dx, dy, sx, sy, snr) = swiftir.swim(remodimg, localimg)
    print('-> ', dx, dy, sx, sy, snr)

S = db.sel('select S from runs where r=%s', (r,))[0][0]

rms = []
for s in range(S):
    rms.append((r,m,s))

db.exe('create table if not exists studystackup0 ( s integer, sd float )')
db.exe('delete from studystackup0')


swiftir.remod(rms, loader=loadtile, saver=swimtile)
