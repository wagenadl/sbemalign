import rawimage
import pyqplot as qp
import numpy as np
import aligndb
import swiftir

db = aligndb.DB()

(r,m1,m2,s,ix1,ix2,iy1,iy2,x1,y1,x2,y2,snrc) \
 = db.sel('''select r,m1,m2,s,ix1,ix2,iy1,iy2,x1,y1,x2,y2,snrc from slicealignq5pos where r=1 and m1=0 and m2=1 and s=3 order by ix1''')[2]

img1 = rawimage.partialq5img(r, m1, s, ix1, iy1)
img2 = rawimage.partialq5img(r, m2, s, ix2, iy2)

qp.figure('/tmp/s1', 8, 4)
qp.subplot(1,2,1)
qp.imsc(img1, xx=np.array([0,684]), yy=[0,-684])
qp.pen('r')
qp.marker('o', 3)
qp.mark(x1, -y1)
qp.shrink()

qp.subplot(1,2,2)
qp.imsc(img2, xx=[0,684], yy=[0,-684])
qp.pen('r')
qp.marker('o', 3)
qp.mark(x2, -y2)
qp.shrink()


