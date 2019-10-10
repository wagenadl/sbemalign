#!/usr/bin/python3

# This sees if the naive shifting made things better

import rawimage
import aligndb
import swiftir

db = aligndb.DB()

r = 21
m = 4
q = 25

def loadtile(rms):
    # Loads a tile with shift from stackup0 table
    r,m,s = rms
    print('loadtile', r, m, s)
    img = rawimage.loadimage(rawimage.scaledtile(r,m,s,q))
    try:
        (dx,dy) = db.sel('select dx,dy from stackup0 where r=%s and m=%s and s=%s',
                     (r,m,s))[0]
    except Exception as e:
        print(e)
        (dx,dy) = (0,0)
    dx /= q
    dy /= q
    (dx,dy)= (0,0)
    Y,X = img.shape
    return swiftir.extractStraightWindow(img, (X/2.-dx, Y/2.-dy), siz=Y)

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
    db.exe('''insert into stackup1 (r,m,s, dx,dy,sx,sy,snr)
    values (%s,%s,%s, %s,%s,%s,%s,%s)''', (r,m,s,
                                           float(dx*q), float(dy*q),
                                           float(sx*q), float(sy*q),
                                           float(snr)))

S = db.sel('select S from runs where r=%s', (r,))[0][0]

db.exe('''create table if not exists stackup1 (
r integer, m integer, s integer,
dx float, dy float,
sx float, sy float,
snr float )''')

db.exe('''delete from stackup1''')


rms = []
for s in range(S):
    rms.append((r,m,s))

swiftir.remod(rms, loader=loadtile, saver=swimtile)
