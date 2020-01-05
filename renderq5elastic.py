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
import config

tbl = 'renderq5elasticdone'
odir = f'{config.sclroot}/q5elastic'

db = aligndb.DB()
ri = db.runinfo()
x0, y0, x1, y1 = renderq5utils.globalbbox()
W = x1
H = y1

def droptable():
    db.nofail(f'drop table {tbl}')

def maketable():
    db.exe(f'''create table if not exists {tbl} (
    z integer
    )''')
    db.exe(f'''create index if not exists {tbl}_z on {tbl} (z)''')

def render(z):
    print(f'Rendering z{z}')
    img = np.zeros((H,W), dtype=np.uint8)
    r, s = ri.findz(z)
    M = ri.nmontages(r)
    for m in range(M):
        print(f'Rendering Z{z} M{m}')
        tile = rawimage.fullq5img(r, m, s)
        if tile is None:
            print(f'Not rendering Z{z} M{m} - No image')
            continue
        xx, yy = renderq5utils.rendergrid(r, m, s)
        IX = len(xx) - 1
        IY = len(yy) - 1
        x0, y0 = renderq5utils.rigidtileposition(r, m, s)
        for ix in range(IX):
            for iy in range(IY):
                #print(f'Rendering Z{z} M{m} IX{ix} IY{iy}')
                xl1 = xx[ix]
                xr1 = xx[ix+1]
                yt1 = yy[iy]
                yb1 = yy[iy+1]
                if xr1 <= xl1 or yb1 <= yt1:
                    continue
                xywh = [xl1, yt1, xr1-xl1, yb1-yt1]
                xmdl = np.array([xx[ix], xx[ix], xx[ix+1], xx[ix+1]])
                ymdl = np.array([yy[iy], yy[iy+1], yy[iy], yy[iy+1]])
                dx, dy = renderq5utils.interpolatedshifts(r, m, s, xmdl, ymdl)
                xtile = xmdl - x0 - dx
                ytile = ymdl - y0 - dy
                warp.warpToRectangle(img, xywh, tile, xmdl, ymdl, xtile, ytile)
    dr = f'{odir}/Z{z//100}'
    if not(os.path.exists(dr)):
        os.mkdir(dr)
    fn = f'{dr}/{z%100}.jpg'
    cv2.imwrite(fn, img)
    db.exe(f'insert into {tbl} (z) values ({z})')
    
def perhapsrender(z):
    if db.sel(f'select count(1) from {tbl} where z={z}')[0][0]==0:
        render(z)

if not(os.path.exists(odir)):
    os.mkdir(odir)
maketable()    
args = sys.argv
args.pop(0)
if len(args)>0:
  for a in args:
    z = int(a)
    render(z)
else: 
  fac = factory.Factory(12)
  for z in range(4,9604):
    fac.request(perhapsrender, z)
  fac.shutdown()
  
