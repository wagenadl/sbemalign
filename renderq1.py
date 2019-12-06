#!/usr/bin/python3

import aligndb
import renderq5utils
import numpy as np
import cv2
import factory
import os
import rawimage
import warp

tbl = 'renderq1done'
odir = '/lsi2/dw/170428/q1pyramid'


db = aligndb.DB()
ri = db.runinfo()
x0, y0, x1, y1 = renderq5utils.globalbbox()
W = x1
H = y1
Q = 5

def droptable():
    db.nofail(f'drop table {tbl}')

def maketable():
    db.exe(f'''create table if not exists {tbl} (
    z integer
    )''')
    db.exe(f'''create index if not exists {tbl}_z on {tbl} (z)''')

def render(z):
    print(f'Rendering z{z}')
    img = np.zeros((H*Q, W*Q), dtype=np.uint8)
    r, s = ri.findz(z)
    M = ri.nmontages(r)
    xx0 = []
    yy0 = []
    for m in range(M):
        print(f'Loading Z{z} M{m}')
        tile = rawimage.fullq1img(r, m, s)
        xx, yy = renderq5utils.rendergrid(r, m, s)
        IX = len(xx) - 1
        IY = len(yy) - 1
        x0, y0 = renderq5utils.rigidtileposition(r, m, s)
        xx0.append(x0)
        yy0.append(y0)
        fac = factory.Factory(7)
        def renderone(ix, iy):
            print(f'Rendering Z{z} M{m} IX{ix} IY{iy}')
            xl1 = xx[ix]
            xr1 = xx[ix+1]
            yt1 = yy[iy]
            yb1 = yy[iy+1]
            xywh = [xl1*Q, yt1*Q, (xr1-xl1)*Q, (yb1-yt1)*Q]
            xmdl = np.array([xx[ix], xx[ix], xx[ix+1], xx[ix+1]])
            ymdl = np.array([yy[iy], yy[iy+1], yy[iy], yy[iy+1]])
            dx, dy = renderq5utils.interpolatedshifts(r, m, s, xmdl, ymdl)
            xtile = xmdl - (x0 + dx)
            ytile = ymdl - (y0 + dy)
            warp.warpToRectangle(img, xywh, tile,
                                 xmdl*Q, ymdl*Q,
                                 xtile*Q, ytile*Q)
        for ix in range(IX):
            for iy in range(IY):
                fac.request(renderone, ix, iy)
        fac.shutdown()
    dr = f'{odir}/Z{z//100}'
    if not os.path.exists(dr):
        os.mkdir(dr)
    dr += f'/{z%100}'
    if not os.path.exists(dr):
        os.mkdir(dr)

    x0 = np.min(xx0)
    y0 = np.min(yy0)
    x1 = np.max(xx0) + 17100
    y1 = np.max(yy0) + 17100
    for a in range(8):
        dr1 = dr + f'/A{a}'
        if not os.path.exists(dr1):
            os.mkdir(dr1)
        x0t = int(np.floor(x0)) // 512
        y0t = int(np.floor(y0)) // 512
        x1t = (int(np.ceil(x1))+511) // 512
        y1t = (int(np.ceil(y1))+511) // 512
        def saveone(x, y): # tile numbers!
            print(f'Saving Z{z} A{a} X{x} Y{y}')
            xl = x*512
            xr = xl + 512
            yt = y*512
            yb = yt + 512
            cv2.imwrite(f'{dr}/A{a}/Y{y}/X{x}.jpg', img[yt:yb, xl:xr])
        fac = factory.Factory(8)
        for yt in range(y0t, y1t):
            dr2 = dr1 + f'/Y{yt}'
            if not os.path.exists(dr2):
                os.mkdir(dr2)
            for xt in range(x0t, x1t):
                fac.request(saveone, xt, yt)
        fac.shutdown()
        x0 = x0 // 2
        y0 = y0 // 2
        x1 = (x1 + 1) // 2
        y1 = (y1 + 1) // 2
        img = rawimage.iscale(img, 2)
    db.exe(f'insert into {tbl} (z) values ({z})')

def perhapsrender(z):
    if db.sel(f'select count(1) from {tbl} where z={z}')[0][0]==0:
        render(z)

if not(os.path.exists(odir)):
    os.mkdir(odir)
maketable()    

fac0 = factory.Factory(4)
for z in range(4,9604):
    fac0.request(perhapsrender, z)
    
