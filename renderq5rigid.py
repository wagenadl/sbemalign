#!/usr/bin/python3

import aligndb
import sys
import os
import numpy as np
import rawimage
import pyqplot as qp
import factory
import swiftir
import cv2

monttbl = 'solveq5mont'
rigidtbl = 'solveq5rigidtile'

X = Y = 684*5 # Full tile size in q5 space!
nz = 200

db = aligndb.DB()
ri = db.runinfo()

def imagespace(z0):
    x0, y0, x1, y1 = db.sel(f'''select 
    min(mo.x+ri.x), min(mo.y+ri.y),
    max(mo.x+ri.x), max(mo.y+ri.y)
    from {monttbl} as mo
    inner join {rigidtbl} as ri on mo.z0=ri.z0 and mo.r=ri.r and mo.m=ri.m
    where mo.z0={z0}''')[0]
    x0 = int(x0)
    y0 = int(y0)
    W = int(x1+X+1)-x0
    H = int(y1+Y+1)-y0
    return (x0,y0,W,H)
    

def renderslice(z0, r, s):
    # Returns an image and x,y coordinates in model space of its top left
    # Positions of tiles
    x0, y0 = db.vsel(f'''select x,y from {monttbl} 
    where z0={z0} and r={r} order by m''')
    x1, y1 = db.vsel(f'''select x,y from {rigidtbl} 
    where z0={z0} and r={r} and s={s} order by m''')

    # Full extent of all tiles
    xxl = x0 + x1
    xxr = xxl + X
    yyt = y0 + y1
    yyb = yyt + Y

    R = ri.nrows(r)
    C = ri.ncolumns(r)

    (x0,y0,W,H) = imagespace(z0)
    img = np.zeros((H,W), dtype=np.uint8) + 255

    if C==1:
        yjoin = np.zeros(R+1, dtype=int)
        yjoin[0] = int(yyt[0]+1)
        yjoin[-1] = int(yyb[-1])
        for row in range(R-1):
            yjoin[row+1] = int((yyb[row] + yyt[row+1])/2 + .5)
        for row in range(R):
            m = row
            print(f'  Loading fullq5img {r} {m} {s}')
            im1 = rawimage.fullq5img(r, m, s)
            if im1 is None:
                print(f'Skipping tile R{r} M{m} S{s} - no image')
                continue
            yt = yjoin[row] # in model coordinates
            yb = yjoin[row+1]
            xl = int(xxl[m]+1)
            xr = xl + X - 2
            yt1 = yjoin[row] - yyt[m] # in tile coordinates
            yb1 = yt1 + yb-yt
            xl1 = xl - xxl[m]
            xr1 = xl1 + xr-xl
            cen = ((xl1+xr1)/2, (yt1+yb1)/2)
            siz = (xr-xl, yb-yt)
            im1 = swiftir.extractStraightWindow(im1, cen, siz)
            img[yt-y0:yb-y0,xl-x0:xr-x0] =  im1
    elif C==2:
        yjoin = np.zeros(R-1, dtype=int)
        xjoin = np.zeros(R, dtype=int)
        for row in range(R-1):
            yjoin[row] = (yyb[row*C] + yyb[row*C+1]
                          + yyt[(row+1)*C] + yyt[(row+1)*C+1]) / 4
        for row in range(R):
            xjoin[row] = (xxr[row*C] + xxl[row*C+1])/2
        for col in range(C):
            for row in range(R):
                m = col + C*row
                print(f'  Loading fullq5img {r} {m} {s}')
                im1 = rawimage.fullq5img(r, m, s)
                if im1 is None:
                    print(f'Skipping tile R{r} M{m} S{s} - no image')
                    continue
 
                if row==0:
                    yt = int(yyt[m]+1)
                else:
                    yt = yjoin[row-1]
                if row==R-1:
                    yb = int(yyb[m])
                else:
                    yb = yjoin[row]
                if col==0:
                    xl = int(xxl[m]+1)
                else:
                    xl = xjoin[row]
                if col==C-1:
                    xr = int(xxr[m])
                else:
                    xr = xjoin[row]
                yt1 = yt - yyt[m]
                yb1 = yt1 + yb-yt
                xl1 = xl - xxl[m]
                xr1 = xl1 + xr-xl
                cen = ((xl1+xr1)/2, (yt1+yb1)/2)
                siz = (xr-xl, yb-yt)
                im1 = swiftir.extractStraightWindow(im1, cen, siz)
                img[yt-y0:yb-y0,xl-x0:xr-x0] = im1
    else:
        raise ValueError('Bad column count')    
    return img, x0, y0

argc = len(sys.argv)
if argc>=2:
    z0 = int(argc[1])
else:
    z0 = None
if argc>=3:
    r = int(argc[2])
else:
    r = None
if argc>=4:
    s = int(argc[3])
else:
    s = None

if s is None:
   fac = factory.Factory(12)

def rendersvslice(z0, r, s):
    odir = f'/lsi2/dw/170428/rigidq5/Z{z0}'
    if not os.path.exists(odir):
        os.mkdir(odir)
    ofn = odir + f'/R{r}S{s}.jpg'
    if not os.path.exists(ofn):
        print(f'Rendering R{r}S{s} in z0={z0}')
        img, x0, y0 = renderslice(z0, r, s)
        rawimage.saveimage(img, ofn)

def rendersvrun(z0, r, s):
    if s is None:
        sv = ri.subvolume(z0, nz)
        for rl in sv:
            if rl.r==r:
                ss = range(rl.s0, rl.s1)
        for s in ss:
            fac.request(rendersvslice, z0, r, s)
    else:
        rendersvslice(z0, r, s)    
            
def rendersubvol(z0, r, s):
    if r is None:
        sv = ri.subvolume(z0, nz)
        for rl in sv:
            rendersvrun(z0, rl.r, None)
    else:
        rendersvrun(z0, r, s)
        
if z0 is None:
    for z0 in range(4,9500,200):
        rendersubvol(z0, None, None)
else:
    rendersubvol(z0, r, s)
       
if s is None:
   fac.shutdown()

