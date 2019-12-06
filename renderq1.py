#!/usr/bin/python3

import aligndb
import renderq5utils
import numpy as np
import cv2
import factory
import os
import rawimage
import warp
import sys

tbl = 'renderq1done'
odir = '/lsi2/dw/170428/q1pyramid'


db = aligndb.DB()
ri = db.runinfo()
x0, y0, x1, y1 = renderq5utils.globalbbox()
W = x1
H = y1
Q = 5
TS = 512 # Tile size
FS = 17100 # Full size

def droptable():
    db.nofail(f'drop table {tbl}')

def maketable():
    db.exe(f'''create table if not exists {tbl} (
    z integer
    )''')
    db.exe(f'''create index if not exists {tbl}_z on {tbl} (z)''')

def render(z):
    print(f'Rendering z{z}')
    X = TS*((W*Q + TS-1) // TS)
    Y = TS*((H*Q + TS-1) // TS)
    img = np.zeros((Y,X), dtype=np.uint8)
    r, s = ri.findz(z)
    M = ri.nmontages(r)
    xx0 = []
    yy0 = []
    for m in range(M):
        print(f'Loading Z{z} M{m}')
        tile = rawimage.fullq1img(r, m, s)
        print(f'Replacing 0 by 1 in Z{z} M{m}')
        tile[tile==0] = 1
        print(f'Processing Z{z} M{m}')
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

    img1 = img
    # Find extent actually filled
    x0 = int(np.floor(np.min(xx0)*Q))
    y0 = int(np.floor(np.min(yy0)*Q))
    x1 = int(np.ceil(np.max(xx0)*Q + FS))
    y1 = int(np.ceil(np.max(yy0)*Q + FS))
    for a in range(9):
        dr1 = dr + f'/A{a}'
        if not os.path.exists(dr1):
            os.mkdir(dr1)
        x0t = x0 // TS
        y0t = y0 // TS
        x1t = (x1+TS-1) // TS
        y1t = (y1+TS-1) // TS
        def saveone(x, y, dr2, img1): # x and y are tile numbers!
            print(f'Saving Z{z} A{a} X{x} Y{y}')
            xl = x*TS
            xr = xl + TS
            yt = y*TS
            yb = yt + TS
            img2 = img1[yt:yb, xl:xr]
            print(f'Z{z} A{a} X{x} Y{y}', xl,yt,img1.shape,
                  img2.shape,np.mean(img2))
            if not cv2.imwrite(f'{dr2}/X{x}.jpg', img2):
                print(f'Failed to save Z{z} A{a} X{x} Y{y}')
                sys.exit(1)
        fac = factory.Factory(8)
        for yt in range(y0t, y1t):
            dr2 = dr1 + f'/Y{yt}'
            if not os.path.exists(dr2):
                os.mkdir(dr2)
            for xt in range(x0t, x1t):
                fac.request(saveone, xt, yt, dr2, img1)
        fac.shutdown()
        x0 = x0 // 2
        y0 = y0 // 2
        x1 = (x1 + 1) // 2
        y1 = (y1 + 1) // 2
        img1 = rawimage.iscale(img1, 2)
        print(f'Z{z} A{a} Image scaled to ', img1.shape)
        img1 = rawimage.ipad(img1, TS)
        print(f'Z{z} A{a} Image padded to ', img1.shape)
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
    
