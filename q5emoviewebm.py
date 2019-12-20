#!/usr/bin/python3

import cv2
import rawimage
import os
import factory

fdir = '/home/wagenaar/q5eframes'
ovfn = '/lsi2/dw/170428/q5emovie.webm'

def makeframe(z, ifn, ofn):
    print(f'Processing Z{z}')
    img = rawimage.loadimage(ifn)
    Y,X = img.shape
    CY = Y//2
    CX = X//2
    W = 1024
    H = 768
    LX = CX - W//2
    TY = CY - H//2
    img = img[TY:TY+H, LX:LX+W]
    cv2.putText(img, f'{z:04}', (10, H-10), cv2.FONT_HERSHEY_SIMPLEX, 1, (0,0,0))
    rawimage.saveimage(img, ofn)
    

def perhapsmakeframe(z, ofn):
    z1 = z//100
    z2 = z%100
    ifn = f'/lsi2/dw/170428/q5elastic/Z{z1}/{z2}.jpg'
    if os.path.exists(ofn):
        if os.path.getmtime(ofn) > os.path.getmtime(ifn):
            return
    makeframe(z, ifn, ofn)
    
if not os.path.exists(fdir):
    os.mkdir(fdir)
if os.path.exists(ovfn):
    os.unlink(ovfn)

fac = factory.Factory(12)
for z0 in range(9600):
    ofn = f'{fdir}/{z0}.jpg'
    z = z0 + 4
    fac.request(perhapsmakeframe, z, ofn)
fac.shutdown()

print('Converting to movie')
os.system(f'ffmpeg -i "{fdir}/%d.jpg"  -c:v libvpx-vp9 -b:v 0 -crf 35 -threads 8 -an {ovfn}')
