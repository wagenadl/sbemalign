#!/usr/bin/python3

# This does a completely naive stacking of a montage with local
# averaged results

import rawimage
import aligndb
import swiftir
import pyqplot as qp

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
    
    Y,X = remodimg.shape
    R = 512
    marg = (X-R)//2
    remodimg = swiftir.extractStraightWindow(remodimg, siz=R)
    remodimg = swiftir.apodize(remodimg)
    localimg = swiftir.extractStraightWindow(localimg, siz=R)
    localimg = swiftir.apodize(localimg)
    (dx, dy, sx, sy, snr) = swiftir.swim(remodimg, localimg)
    print('-> ', dx, dy, sx, sy, snr)
    db.exe('''insert into stackup0 (r,m,s, dx,dy,sx,sy,snr)
    values (%s,%s,%s, %s,%s,%s,%s,%s)''', (r,m,s,
                                           float(dx*q), float(dy*q),
                                           float(sx*q), float(sy*q),
                                           float(snr)))

S = db.sel('select S from runs where r=%s', (r,))[0][0]

db.exe('''create table if not exists stackup0 (
r integer, m integer, s integer,
dx float, dy float,
sx float, sy float,
snr float )''')

db.exe('''delete from stackup0''')


rms = []
for s in range(S):
    rms.append((r,m,s))

swiftir.remod(rms, loader=loadtile, saver=swimtile)
