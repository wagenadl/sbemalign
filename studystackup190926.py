#!/usr/bin/python3

# This does a completely naive stacking of a montage with local
# averaged results

import rawimage
import aligndb
import swiftir

db = aligndb.DB()

r = 21
m = 4
q = 25

def loadtile(r,m,s):
    # Loads a tile without any shift
    return rawimage.loadimage(rawimage.scaledtile(r,m,s,q))

def swimtile(remodimg, r,m,s, localimg):
    # REMODIMG is the leave-one-out average
    # LOCALIMG is the left-out image
    (dx, dy, sx, sy, snr) = SWIM(remodimg, localimg)
    db.exe('''insert into stackup0 (r,m,s, dx,dy,sx,sy,snr)
    values (%s,%s,%s %s,%s,%s,%s,%s)''' % (r,m,s,
                                           float(dx), float(dy),
                                           float(sx), float(sy),
                                           float(snr)))

S = db.sel('select S from runs where r=%s and m=%s' % (r,m))[0][0]

db.exe('''create table if not exists stackup0 (
r integer, m integer, s integer,
dx float, dy float,
sx float, sy float,
snr float )''')

db.exe('''delete from stackup0''')


rms = []
for s in range(s):
    rms.append((r,m,s))

swiftir.remod(rms, loader=loadtile, saver=swimtile)
