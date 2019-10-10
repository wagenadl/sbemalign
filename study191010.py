%aimport rawimage
import pyqplot as qp
import numpy as np

img1=rawimage.partialq5img(21,2,20, 4,2)
img2=rawimage.partialq5img(21,2,20, 4,3)
Y,X = img1.shape

qp.figure('/tmp/s1', 4, 8)
qp.imsc(img1, rect=[0,0,X,Y])
qp.imsc(img2, rect=[0,-Y,X,Y])
qp.shrink()

img1a=rawimage.partialq5img(21,3,20, 0,2)
img2a=rawimage.partialq5img(21,3,20, 0,3)
Y,X = img1.shape

qp.figure('/tmp/s1a', 4, 8)
qp.imsc(img1a, rect=[0,0,X,Y])
qp.imsc(img2a, rect=[0,-Y,X,Y])
qp.shrink()

######################################################################
import swiftir

apo = swiftir.apodize(img1)
apoa = swiftir.apodize(img1a)
(dx, dy, sx, sy, snr) = swiftir.swim(apo, apoa)

win = swiftir.extractStraightWindow(img1, (X/2-dx/2,Y/2-dy/2),(X//2,Y//2))
wina = swiftir.extractStraightWindow(img1a, (X/2+dx/2,Y/2+dy/2),(X//2,Y//2))

qp.figure('/tmp/s2', 8,4)
qp.subplot(1,2,1)
qp.imsc(win)
qp.shrink()
qp.subplot(1,2,2)
qp.imsc(wina)
qp.shrink()

apo_ = swiftir.apodize(win)
apoa_ = swiftir.apodize(wina)
(dx_, dy_, sx_, sy_, snr_) = swiftir.swim(apo_, apoa_)

# Conclusion: pixel (x,y) on the "A" image fits on (x-dx-dx_,y-dy-dy_) of the
# "base" image
# I will create a db table with (r, m1, m2, s, ix1, iy1, ix2, iy2),
# (dx, dy, sx, sy, snr), and (dx_, dy_, sx_, sy_, snr_).
# The logic will be that m2 pixel (x,y) fits to m1 pixel(x-dx-dx_,y-dy-dy_).
# It will be understood that the "_" alignment was centered at
# (X/2-dx/2, Y/2-dy/2) of the "m1" image.
# Although this information is slightly redundant, I think it will be
# convenient to have it all.
# I'll call the table SLICEALIGNQ5.
# I will do the calculation for each of the 5 edge subtiles in Q5.
# Note that I am _not_ scaling the dx, dy, etc. coordinates back up to Q1.
# That's why "Q5" is in the table name.
#
# Completely separately, I will create a table called MONTAGEALIGNQ5 with
# columns (r, m, s, ix, iy), (dx, dy, sx, sy, snr) that positions pixel
# (x,y) of slice S of montage (R, M) at (x-dx, y-dy) of the local average
# of neighboring slices.
# I will do the calculation for the center subtile and each of the 4 corner
# subtiles in Q5. Perhaps also for the 4 "center-of-edge" subtiles, and
# if the thing runs _really_ fast, also for all the other subtiles.
# Note that I am _not_ scaling the dx, dy, etc. coordinates back up to Q1.
