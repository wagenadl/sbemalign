#!/usr/bin/python3

import cv2
import rawimage
import os
import factory
import renderq5utils

gbb = renderq5utils.globalbbox()
Q = 5
fullwidth = gbb[2] * Q
fullheight = gbb[3] * Q
R = 512 # tile size 
X0 = fullwidth//R//2
Y0 = fullheight//R//2

fdir = '/home/wagenaar/q1eframes'
ovfn = '/lsi2/dw/170428/q1emovie.mp4'

def makeframe(z, ifn1, ifn2, ofn):
    print(f'Processing Z{z}')
    z1 = z//100
    z2 = z%100
    img = np.zeros((R*2,R*2), dtype=np.uint8)
    for y in rage(2):
        for x in range(2):
            ifn = f'/lsi2/dw/170428/q5elastic/Z{z1}/{z2}/Y{Y0+y}/{X0+x}.jpg'
            img1 = rawimage.loadimage(ifn)
            Y, X = img.shape
            if Y != R or X != R:
                raise Exception('Mismatching subtile size')
            img[y*R:(y+1)*R, x*R:(x+1)*R] = img1
    Y,X = img.shape
    CY = Y//2
    CX = X//2
    cv2.putText(img, f'{z:04}', (10, Y-10), cv2.FONT_HERSHEY_SIMPLEX, 1, (0,0,0))
    rawimage.saveimage(img, ofn)
        
def perhapsmakeframe(z, ofn):
    z1 = z//100
    z2 = z%100
    ifn = f'/lsi2/dw/170428/q5elastic/Z{z1}/{z2}/Y{Y0}/{X0}.jpg'
    if os.path.exists(ofn):
        if os.path.getmtime(ofn) > os.path.getmtime(ifn):
            return
    makeframe(z, ofn)
    
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
