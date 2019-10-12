import rawimage
import pyqplot as qp
import numpy as np
import aligndb
import swiftir

db = aligndb.DB()

(r,m1,m2,s,ix1,ix2,iy1,iy2,dx0,dy0, \
 dx,dy,dx,sy,snr,dxb,dyb,sxb,syb,snrb,dxc,dyc,sxc,syc,snrc) \
 = db.vsel('select * from slicealignq5 where r=4 and m1=2 and m2=3 and s=0')

