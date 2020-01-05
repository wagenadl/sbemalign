#!/usr/bin/python3

import cv2
import rawimage
import os
import factory
import config

idir = f'{config.sclroot}/q100elastic'
fdir = f'{config.tmproot}/q100eframes'
ovfn = f'{config.sclroot}/q100emovie.mp4'

def makeframe(z, ifn, ofn):
    print(f'Processing Z{z}')
    img = rawimage.loadimage(ifn)
    Y,X = img.shape
    W = X
    H = Y
    cv2.putText(img, f'{z:04}', (10, H-10), 
        cv2.FONT_HERSHEY_SIMPLEX, 1, (0,0,0))
    cv2.putText(img, f'{z:04}', (8, H-12),
        cv2.FONT_HERSHEY_SIMPLEX, 1, (255,255,255))
    rawimage.saveimage(img, ofn)

def perhapsmakeframe(z, ofn):
    z1 = z//100
    z2 = z%100
    ifn = f'{idir}/Z{z1}/{z2}.jpg'
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
os.system(f'ffmpeg -i "{fdir}/%d.jpg"  -c:v libx265 {ovfn}')
