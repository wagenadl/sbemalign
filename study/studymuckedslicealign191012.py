import rawimage
import pyqplot as qp
import numpy as np
import aligndb
import swiftir

db = aligndb.DB()

r=3
m1=3
m2=5
s=1
ix1=4
iy1=4
ix2=4
iy2=0
sidebyside=0
print(f'Working on {r} {m1}:{m2} {s} {ix1},{iy1}:{ix2},{iy2}')
img1 = rawimage.partialq5img(r,m1,s,ix1,iy1)
img2 = rawimage.partialq5img(r,m2,s,ix2,iy2)
Y,X = img1.shape
if sidebyside:
    img1 = img1[:,X//2:]
    img2 = img2[:,:X//2]
    dx0 = X//2
    dy0 = 0
else:
    img1 = img1[Y//2:,:]
    img2 = img2[:Y//2,:]
    dx0 = 0
    dy0 = Y//2
Y,X = img1.shape        
apo1 = swiftir.apodize(img1)
apo2 = swiftir.apodize(img2)
(dx, dy, sx, sy, snr) = swiftir.swim(apo1, apo2)

win1 = swiftir.extractStraightWindow(img1, (X/2-dx/2,Y/2-dy/2),(X//2,Y//2))
win2 = swiftir.extractStraightWindow(img2, (X/2+dx/2,Y/2+dy/2),(X//2,Y//2))
apo1b = swiftir.apodize(win1)
apo2b = swiftir.apodize(win2)
(dxb, dyb, sxb, syb, snrb) = swiftir.swim(apo1b, apo2b)

qp.figure('/tmp/s1', 8,4)
qp.subplot(1,2,1)
qp.imsc(apo1, xx=np.arange(X), yy=np.arange(Y))
qp.pen('r')
qp.marker('o', 3)
qp.mark(X/2 - dx/2, Y/2 + dy/2)
qp.shrink()
qp.subplot(1,2,2)
qp.imsc(apo2, xx=np.arange(X), yy=np.arange(Y))
qp.pen('r')
qp.marker('o', 3)
qp.mark(X/2 + dx/2, Y/2 - dy/2)
qp.shrink()

qp.figure('/tmp/s3', 8,4)
qp.subplot(1,2,1)
qp.imsc(win1)
qp.shrink()
qp.subplot(1,2,2)
qp.imsc(win2)
qp.shrink()


sta = apo1
mov = apo2
wht = -.65
shf = swiftir.alignmentImage(sta, mov, wht)
xy = swiftir.findPeak(shf)
rad = 5
win = swiftir.extractFoldedAroundPoint(shf, xy, rad)
cm = swiftir.centerOfMass(win)

qp.figure('/tmp/sw', 3, 3)
qp.imsc(win)
