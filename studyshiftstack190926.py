#!/usr/bin/python3

import aligndb
import rawimage
db = aligndb.DB()

r = 21
m = 4

x0,y0 = db.vsel('select xc0, yc0 from avgbetapos where r=%s and m=%s' % (r,m))
ss,xx,yy,dx,dy = db.vsel('select s,xc,yc,dx,dy from betapos natural join uctshift where r=%s and m=%s order by s' % (r,m))

q = 25

img = rawimage.q25img(r, m, 100)
