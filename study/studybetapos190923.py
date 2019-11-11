#!/usr/bin/python3


import aligndb
import pyqplot as qp
import numpy as np
import em170428.runinfo
import rawimage

db = aligndb.DB()
#r,m,s,xc,yc = db.vsel('select r,m,s,xc,yc from betapos')

r,m,s,xc,yc = db.vsel('select r,m,s,xc,yc from betapos where r=9 and m=4')

qp.figure('/tmp/s1')
qp.marker(size=2)
qp.mark(s[s<252], xc[s<252])
qp.figure('/tmp/s2')
qp.marker(size=2)
qp.mark(s[s<252], yc[s<252])
print(np.mean(xc[s<252]))
print(np.mean(yc[s<252]))
print(np.std(xc[s<252]))
print(np.std(yc[s<252]))
print(xc[s==252])
print(yc[s==252])

r,m,s,dx,dy = db.vsel('select r,m,s,dx,dy from uctshift where r=9 and m=4')

qp.figure('/tmp/s3')
qp.marker(size=2)
use = s<252
qp.mark(xc[use], dx[use])

qp.figure('/tmp/s4')
qp.marker(size=2)
qp.mark(yc[use], dy[use])

######################################################################

r,m,s, dx,dy,dx1,dy1,dx2,dy2, cx,cy,cx1,cy1,cx2,cy2 \
    = db.vsel('''select
    bp.r,bp.m,bp.s,
    uc.dx,uc.dy,uc1.dx,uc1.dy,uc2.dx,uc2.dy,
    bp.xc,bp.yc,bp1.xc,bp1.yc,bp2.xc,bp2.yc
    from uctshift as uc
    natural join betapos as bp
    inner join betapos as bp1 on bp.r=bp1.r and bp.m=bp1.m and bp.s=bp1.s-1 
    inner join betapos as bp2 on bp.r=bp2.r and bp.m=bp2.m and bp.s=bp2.s+1
    inner join uctshift as uc1 on uc.r=uc1.r and uc.m=uc1.m and uc.s=uc1.s-1
    inner join uctshift as uc2 on uc.r=uc2.r and uc.m=uc2.m and uc.s=uc2.s+1 ''')

dxe = cx - (cx1+cx2)/2
dye = cy - (cy1+cy2)/2
dxu = dx - (dx1+dx2)/2
dyu = dy - (dy1+dy2)/2

qp.figure('/tmp/s5')
qp.marker(size=2)
qp.mark(dxe, dxu)
qp.shrink(1,1)

qp.figure('/tmp/s5')
qp.marker(size=2)
qp.mark(dye, dyu)
qp.shrink(1,1)

dde = np.sqrt(dxe**2 + dye**2)
use = dde>300
np.sum(use)
